#!/usr/bin/env python3
"""Bootstrap 95% confidence intervals for test-set survival C-indices."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import joblib
import numpy as np
import pyarrow.parquet as pq
import torch
from sksurv.exceptions import NoComparablePairException
from sksurv.metrics import concordance_index_censored

from train_explainable_dl import OmicsDeepSurv, make_blocks, predict_risk


def read_frame(path: Path):
    return pq.read_table(path).to_pandas().set_index("sample_key")


def score(event, time_values, risk):
    return float(concordance_index_censored(event.astype(bool), time_values, risk)[0])


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.input_dir.resolve()
    started = time.perf_counter()
    x_test = read_frame(root / "X_test_selected.parquet")
    y_test = read_frame(root / "y_test.parquet").loc[x_test.index]
    cox_bundle = joblib.load(root / "coxph_baseline.joblib")
    cox_values = cox_bundle["scaler"].transform(x_test[cox_bundle["features"]])
    import pandas as pd
    cox_risk = cox_bundle["model"].predict_partial_hazard(
        pd.DataFrame(cox_values, columns=cox_bundle["features"])
    ).to_numpy()
    rsf_bundle = joblib.load(root / "rsf_baseline.joblib")
    rsf_risk = rsf_bundle["model"].predict(
        x_test[rsf_bundle["features"]].to_numpy(dtype=np.float32)
    )
    checkpoint = torch.load(root / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    dl_model = OmicsDeepSurv(checkpoint["block_dims"], checkpoint["dropout"])
    dl_model.load_state_dict(checkpoint["state_dict"])
    dl_values = (
        x_test[checkpoint["ordered_columns"]].to_numpy(dtype=np.float32) - checkpoint["means"]
    ) / checkpoint["scales"]
    dl_risk = predict_risk(dl_model, make_blocks(dl_values, checkpoint["slices"], torch.device("cpu")))
    risks = {"CoxPH": cox_risk, "Random Survival Forest": rsf_risk, "DeepSurv": dl_risk}
    event = y_test["event"].to_numpy()
    time_values = y_test["time"].to_numpy()
    rng = np.random.default_rng(args.random_state)
    distributions = {model: [] for model in risks}
    skipped = 0
    for iteration in range(args.iterations):
        indices = rng.integers(0, len(y_test), size=len(y_test))
        try:
            iteration_scores = {
                model: score(event[indices], time_values[indices], risk[indices])
                for model, risk in risks.items()
            }
        except NoComparablePairException:
            skipped += 1
        else:
            for model, value in iteration_scores.items():
                distributions[model].append(value)
        if (iteration + 1) % 100 == 0:
            print(f"Bootstrap: {iteration + 1}/{args.iterations}")
    results = {}
    for model, risk in risks.items():
        values = np.asarray(distributions[model], dtype=float)
        results[model] = {
            "test_c_index": score(event, time_values, risk),
            "bootstrap_valid_iterations": len(values),
            "ci_95_percentile_lower": float(np.percentile(values, 2.5)),
            "ci_95_percentile_upper": float(np.percentile(values, 97.5)),
            "bootstrap_mean": float(values.mean()),
            "bootstrap_standard_deviation": float(values.std(ddof=1)),
        }
    elapsed = time.perf_counter() - started
    payload = {
        "test_samples": len(y_test),
        "requested_iterations": args.iterations,
        "skipped_no_comparable_pairs": skipped,
        "random_state": args.random_state,
        "models": results,
        "paired_bootstrap_c_indices": {
            model: [float(value) for value in values]
            for model, values in distributions.items()
        },
        "elapsed_seconds": elapsed,
    }
    (root / "bootstrap_ci_results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    lines = [
        "# Bootstrap Test-Set C-index Confidence Intervals",
        "",
        f"Test samples: {len(y_test)}; requested bootstrap iterations: {args.iterations}; "
        f"valid: {args.iterations-skipped}; seed: {args.random_state}.",
        "Percentile interval: 2.5th–97.5th percentiles.",
        "",
        "| Model | Test C-index | Bootstrap mean | 95% CI | Valid resamples |",
        "|---|---:|---:|---:|---:|",
    ]
    for model, values in results.items():
        lines.append(
            f"| {model} | {values['test_c_index']:.6f} | {values['bootstrap_mean']:.6f} | "
            f"[{values['ci_95_percentile_lower']:.6f}, {values['ci_95_percentile_upper']:.6f}] | "
            f"{values['bootstrap_valid_iterations']} |"
        )
    lines.extend(["", f"Elapsed time: {elapsed:.2f} seconds."])
    (root / "bootstrap_ci_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
