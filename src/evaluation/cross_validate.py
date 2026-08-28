#!/usr/bin/env python3
"""Five-fold survival CV with fold-local univariate Cox feature selection."""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from lifelines import CoxPHFitter
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

from select_features import BLOCKS, collapse_ajcc_onehot, rank_block
from train_explainable_dl import (
    OmicsDeepSurv,
    cindex,
    cox_negative_log_likelihood,
    make_blocks,
    predict_risk,
    set_seed,
)


PREFIXES = {"rna": "rna_", "mirna": "mirna_", "methylation": "meth_", "clinical": "clin_"}


def read_frame(path: Path) -> pd.DataFrame:
    return pq.read_table(path).to_pandas().set_index("sample_key")


def train_dl_fold(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_test: pd.DataFrame,
    y_test: pd.DataFrame,
    args,
    fold_seed: int,
) -> tuple[float, int, float]:
    set_seed(fold_seed)
    columns_by_block = {
        block: [c for c in x_train.columns if c.startswith(prefix)]
        for block, prefix in PREFIXES.items()
    }
    ordered = [c for block in PREFIXES for c in columns_by_block[block]]
    slices = {}
    cursor = 0
    for block in PREFIXES:
        end = cursor + len(columns_by_block[block])
        slices[block] = (cursor, end)
        cursor = end
    x_train = x_train[ordered]
    x_test = x_test[ordered]
    fit_idx, val_idx = train_test_split(
        np.arange(len(x_train)),
        test_size=args.validation_size,
        random_state=fold_seed,
        stratify=y_train["event"].to_numpy(),
    )
    fit = x_train.iloc[fit_idx].to_numpy(dtype=np.float32)
    val = x_train.iloc[val_idx].to_numpy(dtype=np.float32)
    test = x_test.to_numpy(dtype=np.float32)
    means = fit.mean(axis=0, dtype=np.float64).astype(np.float32)
    scales = fit.std(axis=0, dtype=np.float64).astype(np.float32)
    scales[scales == 0] = 1.0
    fit, val, test = (fit - means) / scales, (val - means) / scales, (test - means) / scales
    device = torch.device("cpu")
    fit_blocks = make_blocks(fit, slices, device)
    val_blocks = make_blocks(val, slices, device)
    test_blocks = make_blocks(test, slices, device)
    fit_time = torch.tensor(y_train.iloc[fit_idx]["time"].to_numpy(), dtype=torch.float32)
    fit_event = torch.tensor(y_train.iloc[fit_idx]["event"].to_numpy(), dtype=torch.float32)
    model = OmicsDeepSurv(
        {b: len(columns_by_block[b]) for b in PREFIXES}, dropout=args.dl_dropout
    )
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.dl_learning_rate, weight_decay=args.dl_weight_decay
    )
    best_state, best_val, best_epoch, wait = None, -np.inf, 0, 0
    for epoch in range(1, args.dl_epochs + 1):
        model.train()
        optimizer.zero_grad()
        loss = cox_negative_log_likelihood(model(fit_blocks), fit_time, fit_event)
        loss.backward()
        optimizer.step()
        val_risk = predict_risk(model, val_blocks)
        val_score = cindex(
            y_train.iloc[val_idx]["event"].to_numpy(),
            y_train.iloc[val_idx]["time"].to_numpy(),
            val_risk,
        )
        if val_score > best_val + 1e-6:
            best_val, best_epoch, wait = val_score, epoch, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            wait += 1
        if wait >= args.dl_patience:
            break
    model.load_state_dict(best_state)
    score = cindex(
        y_test["event"].to_numpy(), y_test["time"].to_numpy(), predict_risk(model, test_blocks)
    )
    return score, best_epoch, best_val


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--cox-penalizer", type=float, default=0.1)
    parser.add_argument("--cox-max-steps", type=int, default=50)
    parser.add_argument("--rsf-trees", type=int, default=500)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--dl-epochs", type=int, default=300)
    parser.add_argument("--dl-patience", type=int, default=30)
    parser.add_argument("--dl-learning-rate", type=float, default=1e-3)
    parser.add_argument("--dl-weight-decay", type=float, default=1e-4)
    parser.add_argument("--dl-dropout", type=float, default=0.3)
    parser.add_argument("--validation-size", type=float, default=0.15)
    parser.add_argument(
        "--allow-legacy-prefilter",
        action="store_true",
        help="Run the historical fixed-10k representation despite the documented outer-CV leakage.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.allow_legacy_prefilter:
        raise RuntimeError(
            "Bu eski script X_train/X_test içindeki sabit methylation-10k filtresini kullanır ve "
            "outer-CV leakage içerir. Düzeltilmiş fold-local analiz için select_features_v2.py "
            "--counts 200 çalıştırın. Eski sonucu yalnız yeniden üretmek için "
            "--allow-legacy-prefilter verilebilir."
        )
    root = args.input_dir.resolve()
    x_train_original = read_frame(root / "X_train.parquet")
    x_test_original = read_frame(root / "X_test.parquet")
    y_train_original = read_frame(root / "y_train.parquet")
    y_test_original = read_frame(root / "y_test.parquet")
    x = pd.concat([x_train_original, x_test_original], axis=0)
    y = pd.concat([y_train_original, y_test_original], axis=0).loc[x.index]
    x = collapse_ajcc_onehot(x, print_categories=True)
    if not x.index.is_unique:
        raise ValueError("Tam CV verisinde yinelenen sample_key var")
    print(f"Tam CV verisi: X={x.shape}, y={y.shape}, event={int(y['event'].sum())}")
    splitter = StratifiedKFold(
        n_splits=args.folds, shuffle=True, random_state=args.random_state
    )
    fold_results = []
    overall_started = time.perf_counter()
    for fold, (train_idx, test_idx) in enumerate(
        splitter.split(np.zeros(len(y)), y["event"].to_numpy()), 1
    ):
        fold_started = time.perf_counter()
        fold_train_x, fold_test_x = x.iloc[train_idx], x.iloc[test_idx]
        fold_train_y, fold_test_y = y.iloc[train_idx], y.iloc[test_idx]
        print(
            f"\nFold {fold}/{args.folds}: train={len(train_idx)}, test={len(test_idx)}, "
            f"train event={int(fold_train_y['event'].sum())}, test event={int(fold_test_y['event'].sum())}"
        )
        selection_started = time.perf_counter()
        selected = []
        block_stats = {}
        for block, prefix in BLOCKS.items():
            candidates = [c for c in x.columns if c.startswith(prefix)]
            ranking = rank_block(
                fold_train_x,
                candidates,
                fold_train_y,
                args.chunk_size,
                50,
                1e-7,
            )
            valid = ranking.dropna(subset=["p_value"])
            chosen = valid.head(min(args.top_n, len(valid)))["feature"].tolist()
            selected.extend(chosen)
            block_stats[block] = {
                "candidates": len(candidates),
                "valid": len(valid),
                "selected": len(chosen),
                "boundary_p_value": float(valid.iloc[len(chosen) - 1]["p_value"]),
            }
        clinical = [c for c in x.columns if c.startswith("clin_")]
        selected.extend(clinical)
        train_selected = fold_train_x[selected]
        test_selected = fold_test_x[selected]
        selection_seconds = time.perf_counter() - selection_started
        print(f"Fold {fold}: seçim tamamlandı; feature={len(selected)}, süre={selection_seconds:.2f} sn")

        event_mask = fold_train_y["event"].astype(bool)
        event_var = train_selected.loc[event_mask].var(axis=0)
        censored_var = train_selected.loc[~event_mask].var(axis=0)
        excluded = train_selected.columns[(event_var <= 1e-12) | (censored_var <= 1e-12)].tolist()
        cox_features = [c for c in selected if c not in excluded]
        scaler = StandardScaler()
        train_scaled = scaler.fit_transform(train_selected[cox_features])
        test_scaled = scaler.transform(test_selected[cox_features])
        cox_data = pd.DataFrame(train_scaled, columns=cox_features)
        cox_data["time"] = fold_train_y["time"].to_numpy()
        cox_data["event"] = fold_train_y["event"].to_numpy()
        cox_started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cox = CoxPHFitter(penalizer=args.cox_penalizer)
            cox.fit(
                cox_data,
                duration_col="time",
                event_col="event",
                fit_options={"max_steps": args.cox_max_steps},
            )
        cox_risk = cox.predict_partial_hazard(
            pd.DataFrame(test_scaled, columns=cox_features)
        ).to_numpy()
        cox_score = float(
            concordance_index_censored(
                fold_test_y["event"].astype(bool).to_numpy(),
                fold_test_y["time"].to_numpy(),
                cox_risk,
            )[0]
        )
        cox_seconds = time.perf_counter() - cox_started

        rsf_started = time.perf_counter()
        rsf = RandomSurvivalForest(
            n_estimators=args.rsf_trees,
            min_samples_split=10,
            min_samples_leaf=5,
            random_state=args.random_state + fold,
            n_jobs=-1,
        )
        rsf.fit(
            train_selected.to_numpy(dtype=np.float32),
            Surv.from_arrays(
                fold_train_y["event"].astype(bool).to_numpy(),
                fold_train_y["time"].to_numpy(),
            ),
        )
        rsf_score = float(
            concordance_index_censored(
                fold_test_y["event"].astype(bool).to_numpy(),
                fold_test_y["time"].to_numpy(),
                rsf.predict(test_selected.to_numpy(dtype=np.float32)),
            )[0]
        )
        rsf_seconds = time.perf_counter() - rsf_started

        dl_started = time.perf_counter()
        dl_score, dl_epoch, dl_val = train_dl_fold(
            train_selected,
            fold_train_y,
            test_selected,
            fold_test_y,
            args,
            args.random_state + fold,
        )
        dl_seconds = time.perf_counter() - dl_started
        fold_seconds = time.perf_counter() - fold_started
        result = {
            "fold": fold,
            "train_samples": len(train_idx),
            "test_samples": len(test_idx),
            "train_events": int(fold_train_y["event"].sum()),
            "test_events": int(fold_test_y["event"].sum()),
            "selected_features": len(selected),
            "selection": block_stats,
            "coxph": {
                "c_index": cox_score,
                "seconds": cox_seconds,
                "warnings": [str(w.message) for w in caught],
                "excluded_separation_features": excluded,
            },
            "rsf": {"c_index": rsf_score, "seconds": rsf_seconds},
            "deepsurv": {
                "c_index": dl_score,
                "seconds": dl_seconds,
                "best_epoch": dl_epoch,
                "best_validation_c_index": dl_val,
            },
            "selection_seconds": selection_seconds,
            "total_seconds": fold_seconds,
        }
        fold_results.append(result)
        print(
            f"Fold {fold} sonuç: Cox={cox_score:.6f}, RSF={rsf_score:.6f}, "
            f"DeepSurv={dl_score:.6f}; toplam={fold_seconds:.2f} sn"
        )

    summary = {}
    for model in ("coxph", "rsf", "deepsurv"):
        scores = np.array([f[model]["c_index"] for f in fold_results])
        summary[model] = {
            "mean_c_index": float(scores.mean()),
            "sample_standard_deviation": float(scores.std(ddof=1)),
            "fold_c_indices": scores.tolist(),
        }
    elapsed = time.perf_counter() - overall_started
    payload = {
        "n_samples": len(x),
        "folds": args.folds,
        "stratified_by": "event",
        "random_state": args.random_state,
        "top_n_per_omics_per_fold": args.top_n,
        "fold_results": fold_results,
        "summary": summary,
        "total_seconds": elapsed,
    }
    (root / "cv_results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    lines = [
        "# Five-Fold Stratified Cross-Validation Results",
        "",
        f"Samples: {len(x)}; folds: {args.folds}; stratification: event; seed: {args.random_state}.",
        "Feature selection was recomputed using only each outer fold's training subset.",
        "",
        "| Fold | Train/Test | CoxPH C-index | RSF C-index | DeepSurv C-index | Selection s | Total s |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for f in fold_results:
        lines.append(
            f"| {f['fold']} | {f['train_samples']}/{f['test_samples']} | "
            f"{f['coxph']['c_index']:.6f} | {f['rsf']['c_index']:.6f} | "
            f"{f['deepsurv']['c_index']:.6f} | {f['selection_seconds']:.2f} | {f['total_seconds']:.2f} |"
        )
    lines.extend(["", "## Mean ± sample standard deviation", ""])
    for model, label in (("coxph", "CoxPH"), ("rsf", "Random Survival Forest"), ("deepsurv", "DeepSurv")):
        lines.append(
            f"- {label}: {summary[model]['mean_c_index']:.6f} ± "
            f"{summary[model]['sample_standard_deviation']:.6f}"
        )
    lines.extend(["", f"Total elapsed time: {elapsed:.2f} seconds."])
    (root / "cv_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
