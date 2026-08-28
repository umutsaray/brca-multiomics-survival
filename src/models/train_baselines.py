#!/usr/bin/env python3
"""Train multivariable CoxPH and Random Survival Forest baselines."""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
import pyarrow.parquet as pq
from lifelines import CoxPHFitter
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from sklearn.preprocessing import StandardScaler


def read_frame(path: Path):
    return pq.read_table(path).to_pandas().set_index("sample_key")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--cox-penalizer", type=float, default=0.1)
    parser.add_argument("--cox-max-steps", type=int, default=50)
    parser.add_argument("--rsf-trees", type=int, default=500)
    parser.add_argument("--rsf-min-split", type=int, default=10)
    parser.add_argument("--rsf-min-leaf", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=-1)
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.input_dir.resolve()
    previous_results = None
    previous_path = root / "baseline_results.json"
    if previous_path.is_file():
        previous_results = json.loads(previous_path.read_text(encoding="utf-8"))
    x_train = read_frame(root / "X_train_selected.parquet")
    x_test = read_frame(root / "X_test_selected.parquet")
    y_train = read_frame(root / "y_train.parquet").loc[x_train.index]
    y_test = read_frame(root / "y_test.parquet").loc[x_test.index]
    print(f"Baseline girdileri: X_train={x_train.shape}, X_test={x_test.shape}")
    event_mask = y_train["event"].astype(bool)
    conditional_variances_event = x_train.loc[event_mask].var(axis=0)
    conditional_variances_censored = x_train.loc[~event_mask].var(axis=0)
    separation_columns = x_train.columns[
        (conditional_variances_event <= 1e-12) | (conditional_variances_censored <= 1e-12)
    ].tolist()
    cox_features = [c for c in x_train.columns if c not in separation_columns]
    print(
        f"CoxPH complete-separation ön kontrolü: {len(separation_columns)} sütun dışlandı; "
        f"kalan={len(cox_features)}"
    )
    if separation_columns:
        print(f"CoxPH dışlanan sütunlar: {separation_columns}")
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(x_train[cox_features]).astype("float64")
    test_scaled = scaler.transform(x_test[cox_features]).astype("float64")

    cox_started = time.perf_counter()
    cox_data = __import__("pandas").DataFrame(train_scaled, columns=cox_features)
    cox_data["time"] = y_train["time"].to_numpy()
    cox_data["event"] = y_train["event"].to_numpy()
    cox = CoxPHFitter(penalizer=args.cox_penalizer)
    with warnings.catch_warnings(record=True) as caught_cox_warnings:
        warnings.simplefilter("always")
        cox.fit(
            cox_data,
            duration_col="time",
            event_col="event",
            show_progress=False,
            fit_options={"max_steps": args.cox_max_steps},
        )
    cox_warning_messages = [str(item.message) for item in caught_cox_warnings]
    for message in cox_warning_messages:
        print(f"UYARI CoxPH: {message}")
    cox_risk = cox.predict_partial_hazard(
        __import__("pandas").DataFrame(test_scaled, columns=cox_features)
    ).to_numpy()
    cox_cindex = float(
        concordance_index_censored(
            y_test["event"].astype(bool).to_numpy(), y_test["time"].to_numpy(), cox_risk
        )[0]
    )
    cox_seconds = time.perf_counter() - cox_started
    joblib.dump({"model": cox, "scaler": scaler, "features": cox_features}, root / "coxph_baseline.joblib")
    print(f"CoxPH test C-index: {cox_cindex:.6f}; süre: {cox_seconds:.2f} sn")

    survival_train = Surv.from_arrays(
        event=y_train["event"].astype(bool).to_numpy(), time=y_train["time"].to_numpy()
    )
    rsf_started = time.perf_counter()
    rsf = RandomSurvivalForest(
        n_estimators=args.rsf_trees,
        min_samples_split=args.rsf_min_split,
        min_samples_leaf=args.rsf_min_leaf,
        random_state=args.random_state,
        n_jobs=args.n_jobs,
    )
    rsf.fit(x_train.to_numpy(dtype=np.float32), survival_train)
    rsf_risk = rsf.predict(x_test.to_numpy(dtype=np.float32))
    rsf_cindex = float(
        concordance_index_censored(
            y_test["event"].astype(bool).to_numpy(), y_test["time"].to_numpy(), rsf_risk
        )[0]
    )
    rsf_seconds = time.perf_counter() - rsf_started
    joblib.dump({"model": rsf, "features": x_train.columns.tolist()}, root / "rsf_baseline.joblib")
    print(f"RSF test C-index: {rsf_cindex:.6f}; süre: {rsf_seconds:.2f} sn")

    results = {
        "train_shape": list(x_train.shape),
        "test_shape": list(x_test.shape),
        "models": {
            "coxph": {
                "test_c_index": cox_cindex,
                "penalizer": args.cox_penalizer,
                "max_steps": args.cox_max_steps,
                "warnings": cox_warning_messages,
                "complete_separation_columns_excluded": separation_columns,
                "feature_count_after_exclusion": len(cox_features),
                "elapsed_seconds": cox_seconds,
                "model_file": "coxph_baseline.joblib",
            },
            "random_survival_forest": {
                "test_c_index": rsf_cindex,
                "n_estimators": args.rsf_trees,
                "min_samples_split": args.rsf_min_split,
                "min_samples_leaf": args.rsf_min_leaf,
                "random_state": args.random_state,
                "elapsed_seconds": rsf_seconds,
                "model_file": "rsf_baseline.joblib",
            },
        },
    }
    if previous_results and "models" in previous_results:
        prior_comparison = previous_results.get("comparison_to_previous_run", {})
        previous_cox = prior_comparison.get(
            "previous_coxph_test_c_index",
            previous_results["models"]["coxph"]["test_c_index"],
        )
        previous_rsf = prior_comparison.get(
            "previous_rsf_test_c_index",
            previous_results["models"]["random_survival_forest"]["test_c_index"],
        )
        results["comparison_to_previous_run"] = {
            "previous_coxph_test_c_index": previous_cox,
            "new_coxph_test_c_index": cox_cindex,
            "coxph_delta": cox_cindex - previous_cox,
            "previous_rsf_test_c_index": previous_rsf,
            "new_rsf_test_c_index": rsf_cindex,
            "rsf_delta": rsf_cindex - previous_rsf,
        }
        print(
            f"Önceki -> yeni CoxPH C-index: {previous_cox:.6f} -> {cox_cindex:.6f} "
            f"(delta={cox_cindex-previous_cox:+.6f})"
        )
        print(
            f"Önceki -> yeni RSF C-index: {previous_rsf:.6f} -> {rsf_cindex:.6f} "
            f"(delta={rsf_cindex-previous_rsf:+.6f})"
        )
    (root / "baseline_results.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
