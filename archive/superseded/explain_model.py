#!/usr/bin/env python3
"""Explain the trained omics-block DeepSurv risk score with SHAP gradients."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pyarrow.parquet as pq
import shap
import torch

from train_explainable_dl import FlatOmicsWrapper, OmicsDeepSurv


COLORS = {
    "rna": "#2563EB",
    "mirna": "#16A34A",
    "methylation": "#EA580C",
    "clinical": "#7C3AED",
}


def read_frame(path: Path):
    return pq.read_table(path).to_pandas().set_index("sample_key")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--explain-samples", type=int, default=75)
    parser.add_argument("--background-samples", type=int, default=50)
    parser.add_argument("--random-state", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    np.random.seed(args.random_state)
    torch.manual_seed(args.random_state)
    root = args.input_dir.resolve()
    started = time.perf_counter()
    checkpoint = torch.load(root / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    model = OmicsDeepSurv(checkpoint["block_dims"], checkpoint["dropout"])
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    wrapper = FlatOmicsWrapper(model, checkpoint["slices"])
    wrapper.eval()
    x_train = read_frame(root / "X_train_selected.parquet")[checkpoint["ordered_columns"]]
    x_test = read_frame(root / "X_test_selected.parquet")[checkpoint["ordered_columns"]]
    train_values = (
        x_train.to_numpy(dtype=np.float32) - checkpoint["means"]
    ) / checkpoint["scales"]
    test_values = (
        x_test.to_numpy(dtype=np.float32) - checkpoint["means"]
    ) / checkpoint["scales"]
    rng = np.random.default_rng(args.random_state)
    background_n = min(args.background_samples, len(train_values))
    explain_n = min(args.explain_samples, len(test_values))
    background_idx = rng.choice(len(train_values), size=background_n, replace=False)
    explain_idx = rng.choice(len(test_values), size=explain_n, replace=False)
    background = torch.from_numpy(train_values[background_idx])
    explained = torch.from_numpy(test_values[explain_idx])
    print(
        f"SHAP: background={background_n}, açıklanan test örneği={explain_n}, "
        f"özellik={test_values.shape[1]}"
    )
    explainer = shap.GradientExplainer(wrapper, background)
    shap_values = explainer.shap_values(explained)
    if isinstance(shap_values, list):
        shap_values = shap_values[0]
    shap_values = np.asarray(shap_values)
    if shap_values.ndim == 3:
        shap_values = shap_values[..., 0]
    importance = np.mean(np.abs(shap_values), axis=0)
    features = checkpoint["ordered_columns"]
    block_for_feature = []
    for feature in features:
        if feature.startswith("rna_"):
            block_for_feature.append("rna")
        elif feature.startswith("mirna_"):
            block_for_feature.append("mirna")
        elif feature.startswith("meth_"):
            block_for_feature.append("methylation")
        else:
            block_for_feature.append("clinical")
    block_totals = {
        block: float(importance[np.array(block_for_feature) == block].sum())
        for block in COLORS
    }
    top_indices = np.argsort(importance)[::-1][:20]
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    blocks = list(COLORS)
    axes[0].bar(blocks, [block_totals[b] for b in blocks], color=[COLORS[b] for b in blocks])
    axes[0].set_title("Total mean |SHAP| by omics block")
    axes[0].set_ylabel("Sum of mean absolute SHAP values")
    axes[0].tick_params(axis="x", rotation=25)
    top_names = [features[i] for i in top_indices][::-1]
    top_values = importance[top_indices][::-1]
    top_blocks = [block_for_feature[i] for i in top_indices][::-1]
    axes[1].barh(top_names, top_values, color=[COLORS[b] for b in top_blocks])
    axes[1].set_title("Top 20 features by mean |SHAP|")
    axes[1].set_xlabel("Mean absolute SHAP value")
    axes[1].legend(
        handles=[Patch(facecolor=color, label=block) for block, color in COLORS.items()],
        loc="lower right",
    )
    fig.tight_layout()
    output_png = root / "shap_omics_summary.png"
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    lines = [
        "# SHAP Top Features",
        "",
        f"Explained test samples: {explain_n}; background train samples: {background_n}.",
        "",
        "| Rank | Feature | Omics block | Mean absolute SHAP |",
        "|---:|---|---|---:|",
    ]
    for rank, idx in enumerate(top_indices, 1):
        lines.append(
            f"| {rank} | `{features[idx]}` | {block_for_feature[idx]} | {importance[idx]:.10g} |"
        )
    lines.extend(["", "## Block-level total importance", ""])
    for block in blocks:
        lines.append(f"- {block}: {block_totals[block]:.10g}")
    elapsed = time.perf_counter() - started
    lines.extend(["", f"SHAP elapsed time: {elapsed:.2f} seconds."])
    (root / "shap_top_features.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"SHAP tamamlandı: {elapsed:.2f} sn; grafik={output_png.name}")
    print("İlk 20 özellik:")
    for rank, idx in enumerate(top_indices, 1):
        print(f"  {rank}. {features[idx]} ({block_for_feature[idx]}): {importance[idx]:.10g}")


if __name__ == "__main__":
    main()
