#!/usr/bin/env python3
"""Train an omics-block DeepSurv-style PyTorch model."""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import torch
from sklearn.model_selection import train_test_split
from sksurv.metrics import concordance_index_censored
from torch import nn


PREFIXES = {"rna": "rna_", "mirna": "mirna_", "methylation": "meth_", "clinical": "clin_"}


def read_frame(path: Path):
    return pq.read_table(path).to_pandas().set_index("sample_key")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


class Encoder(nn.Module):
    def __init__(self, input_dim: int, hidden1: int, hidden2: int, dropout: float):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden1),
            nn.BatchNorm1d(hidden1),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden1, hidden2),
            nn.BatchNorm1d(hidden2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

    def forward(self, x):
        return self.network(x)


class OmicsDeepSurv(nn.Module):
    def __init__(self, block_dims: dict[str, int], dropout: float = 0.3):
        super().__init__()
        hidden = {
            "rna": (128, 32),
            "mirna": (64, 16),
            "methylation": (128, 32),
            "clinical": (32, 8),
        }
        self.block_order = list(PREFIXES)
        self.encoders = nn.ModuleDict(
            {
                block: Encoder(block_dims[block], hidden[block][0], hidden[block][1], dropout)
                for block in self.block_order
            }
        )
        total = sum(hidden[b][1] for b in self.block_order)
        self.head = nn.Sequential(
            nn.Linear(total, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
        )

    def forward(self, blocks: dict[str, torch.Tensor]):
        encoded = [self.encoders[b](blocks[b]) for b in self.block_order]
        return self.head(torch.cat(encoded, dim=1)).squeeze(1)


class FlatOmicsWrapper(nn.Module):
    """Single-tensor wrapper used by SHAP."""
    def __init__(self, model: OmicsDeepSurv, slices: dict[str, tuple[int, int]]):
        super().__init__()
        self.model = model
        self.slices = slices

    def forward(self, x):
        return self.model({b: x[:, start:end] for b, (start, end) in self.slices.items()}).unsqueeze(1)


def cox_negative_log_likelihood(risk, time_values, event_values):
    """Breslow negative partial log-likelihood over the current batch risk set."""
    order = torch.argsort(time_values, descending=True, stable=True)
    risk = risk[order]
    time_values = time_values[order]
    event_values = event_values[order]
    log_cumulative_risk = torch.logcumsumexp(risk, dim=0)
    # For tied times, every event must use the full tied risk set (last row in tie).
    group_end = torch.arange(len(time_values), device=time_values.device)
    for i in range(len(time_values) - 2, -1, -1):
        if time_values[i] == time_values[i + 1]:
            group_end[i] = group_end[i + 1]
    contributions = risk - log_cumulative_risk[group_end]
    events = event_values > 0
    return -contributions[events].mean()


def make_blocks(array: np.ndarray, slices: dict[str, tuple[int, int]], device):
    tensor = torch.as_tensor(array, dtype=torch.float32, device=device)
    return {b: tensor[:, start:end] for b, (start, end) in slices.items()}


def predict_risk(model, blocks):
    model.eval()
    with torch.no_grad():
        return model(blocks).cpu().numpy()


def cindex(event, time_values, risk):
    return float(concordance_index_censored(event.astype(bool), time_values, risk)[0])


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=30)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--validation-size", type=float, default=0.15)
    parser.add_argument("--random-state", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.random_state)
    root = args.input_dir.resolve()
    started = time.perf_counter()
    x_train = read_frame(root / "X_train_selected.parquet")
    x_test = read_frame(root / "X_test_selected.parquet")
    y_train = read_frame(root / "y_train.parquet").loc[x_train.index]
    y_test = read_frame(root / "y_test.parquet").loc[x_test.index]
    columns_by_block = {
        block: [c for c in x_train.columns if c.startswith(prefix)]
        for block, prefix in PREFIXES.items()
    }
    ordered_columns = [c for b in PREFIXES for c in columns_by_block[b]]
    if set(ordered_columns) != set(x_train.columns):
        raise ValueError("Tanınmayan feature prefix'i bulundu")
    x_train = x_train[ordered_columns]
    x_test = x_test[ordered_columns]
    slices = {}
    cursor = 0
    for block in PREFIXES:
        end = cursor + len(columns_by_block[block])
        slices[block] = (cursor, end)
        cursor = end
    print(f"DL girdileri: X_train={x_train.shape}, X_test={x_test.shape}")
    print(f"Blok boyutları: { {b: len(c) for b, c in columns_by_block.items()} }")

    indices = np.arange(len(x_train))
    fit_idx, val_idx = train_test_split(
        indices,
        test_size=args.validation_size,
        random_state=args.random_state,
        stratify=y_train["event"].to_numpy(),
    )
    fit_values = x_train.iloc[fit_idx].to_numpy(dtype=np.float32)
    val_values = x_train.iloc[val_idx].to_numpy(dtype=np.float32)
    test_values = x_test.to_numpy(dtype=np.float32)
    means = fit_values.mean(axis=0, dtype=np.float64).astype(np.float32)
    scales = fit_values.std(axis=0, dtype=np.float64).astype(np.float32)
    scales[scales == 0] = 1.0
    fit_values = (fit_values - means) / scales
    val_values = (val_values - means) / scales
    test_values = (test_values - means) / scales
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Cihaz: {device}; fit={len(fit_idx)}, validation={len(val_idx)}, test={len(x_test)}")
    fit_blocks = make_blocks(fit_values, slices, device)
    val_blocks = make_blocks(val_values, slices, device)
    test_blocks = make_blocks(test_values, slices, device)
    fit_time = torch.as_tensor(y_train.iloc[fit_idx]["time"].to_numpy(), dtype=torch.float32, device=device)
    fit_event = torch.as_tensor(y_train.iloc[fit_idx]["event"].to_numpy(), dtype=torch.float32, device=device)
    val_event = y_train.iloc[val_idx]["event"].to_numpy()
    val_time = y_train.iloc[val_idx]["time"].to_numpy()

    model = OmicsDeepSurv(
        {b: len(columns_by_block[b]) for b in PREFIXES}, dropout=args.dropout
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay
    )
    best_state = None
    best_val = -np.inf
    best_epoch = 0
    wait = 0
    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad()
        risk = model(fit_blocks)
        loss = cox_negative_log_likelihood(risk, fit_time, fit_event)
        loss.backward()
        optimizer.step()
        val_risk = predict_risk(model, val_blocks)
        val_cindex = cindex(val_event, val_time, val_risk)
        if val_cindex > best_val + 1e-6:
            best_val = val_cindex
            best_epoch = epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            wait = 0
        else:
            wait += 1
        if epoch == 1 or epoch % 10 == 0:
            print(
                f"epoch={epoch} loss={loss.item():.6f} val_cindex={val_cindex:.6f} "
                f"best={best_val:.6f}"
            )
        if wait >= args.patience:
            print(f"Early stopping: epoch={epoch}; best epoch={best_epoch}")
            break
    model.load_state_dict(best_state)
    test_risk = predict_risk(model, test_blocks)
    test_cindex = cindex(y_test["event"].to_numpy(), y_test["time"].to_numpy(), test_risk)
    elapsed = time.perf_counter() - started
    print(f"DL test C-index: {test_cindex:.6f}; toplam süre: {elapsed:.2f} sn")

    checkpoint = {
        "state_dict": best_state,
        "block_dims": {b: len(columns_by_block[b]) for b in PREFIXES},
        "columns_by_block": columns_by_block,
        "ordered_columns": ordered_columns,
        "slices": slices,
        "means": means,
        "scales": scales,
        "dropout": args.dropout,
        "best_epoch": best_epoch,
        "best_validation_c_index": best_val,
        "test_c_index": test_cindex,
        "random_state": args.random_state,
    }
    torch.save(checkpoint, root / "explainable_dl_model.pt")
    baseline = json.loads((root / "baseline_results.json").read_text(encoding="utf-8"))
    cox_score = baseline["models"]["coxph"]["test_c_index"]
    rsf_score = baseline["models"]["random_survival_forest"]["test_c_index"]
    comparison = [
        "# Model Comparison",
        "",
        "| Model | Test C-index |",
        "|---|---:|",
        f"| Penalized multivariable CoxPH | {cox_score:.6f} |",
        f"| Random Survival Forest | {rsf_score:.6f} |",
        f"| Omics-block DeepSurv | {test_cindex:.6f} |",
        "",
        f"DL best validation C-index: {best_val:.6f} (epoch {best_epoch}).",
        f"DL elapsed time: {elapsed:.2f} seconds.",
    ]
    (root / "model_comparison.md").write_text("\n".join(comparison) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
