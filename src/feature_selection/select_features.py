#!/usr/bin/env python3
"""Leakage-safe univariate CoxPH feature selection by omics block."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.stats import norm


BLOCKS = {"rna": "rna_", "mirna": "mirna_", "methylation": "meth_"}
AJCC_PREFIX = "clin_categorical__ajcc_pathologic_stage.diagnoses_"
AJCC_GROUPS = {
    "Stage I": ["Stage I", "Stage IA", "Stage IB"],
    "Stage II": ["Stage II", "Stage IIA", "Stage IIB"],
    "Stage III": ["Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"],
    "Stage IV": ["Stage IV"],
}


def read_frame(path: Path) -> pd.DataFrame:
    return pq.read_table(path).to_pandas().set_index("sample_key")


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    table = pa.Table.from_pandas(frame.reset_index(), preserve_index=False).replace_schema_metadata(None)
    pq.write_table(table, path, compression="snappy")


def collapse_ajcc_onehot(frame: pd.DataFrame, *, print_categories: bool = True) -> pd.DataFrame:
    """Collapse verified AJCC one-hot levels; retain unmapped real levels."""
    stage_columns = [c for c in frame.columns if c.startswith(AJCC_PREFIX)]
    actual_categories = [c.removeprefix(AJCC_PREFIX) for c in stage_columns]
    if print_categories:
        print(f"Dosyada görülen gerçek AJCC kategorileri: {actual_categories}")
    expected = {level for levels in AJCC_GROUPS.values() for level in levels}
    missing = sorted(expected.difference(actual_categories))
    if missing:
        raise ValueError(f"AJCC birleştirmesi için beklenen gerçek kategoriler yok: {missing}")
    result = frame.drop(columns=stage_columns).copy()
    for group, levels in AJCC_GROUPS.items():
        source = [f"{AJCC_PREFIX}{level}" for level in levels]
        result[f"{AJCC_PREFIX}{group}"] = frame[source].max(axis=1).astype("float32")
    unmapped = sorted(set(actual_categories).difference(expected))
    for level in unmapped:
        result[f"{AJCC_PREFIX}{level}"] = frame[f"{AJCC_PREFIX}{level}"].astype("float32")
    if print_categories:
        print(
            f"AJCC birleşimi: {len(stage_columns)} dummy -> "
            f"{len(AJCC_GROUPS) + len(unmapped)} dummy; ayrı tutulan eşleme dışı: {unmapped}"
        )
    return result


def univariate_cox_chunk(
    values: np.ndarray,
    time_values: np.ndarray,
    event_values: np.ndarray,
    max_iter: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit independent one-covariate CoxPH models using Breslow ties."""
    means = values.mean(axis=0, dtype=np.float64)
    scales = values.std(axis=0, dtype=np.float64)
    valid = np.isfinite(scales) & (scales > 0)
    x = np.zeros(values.shape, dtype=np.float64)
    x[:, valid] = (values[:, valid] - means[valid]) / scales[valid]
    order = np.argsort(-time_values, kind="stable")
    times_sorted = time_values[order]
    events = event_values[order].astype(bool)
    x = x[order]
    group_end = np.empty(len(times_sorted), dtype=np.int64)
    start = 0
    while start < len(times_sorted):
        end = start
        while end + 1 < len(times_sorted) and times_sorted[end + 1] == times_sorted[start]:
            end += 1
        group_end[start : end + 1] = end
        start = end + 1
    event_rows = np.flatnonzero(events)
    risk_rows = group_end[event_rows]
    event_x_sum = x[event_rows].sum(axis=0)
    beta = np.zeros(values.shape[1], dtype=np.float64)
    information = np.full(values.shape[1], np.nan)
    for _ in range(max_iter):
        eta = np.clip(x * beta, -50.0, 50.0)
        weights = np.exp(eta)
        s0 = np.cumsum(weights, axis=0)[risk_rows]
        s1 = np.cumsum(weights * x, axis=0)[risk_rows]
        s2 = np.cumsum(weights * x * x, axis=0)[risk_rows]
        score = event_x_sum - (s1 / s0).sum(axis=0)
        information = (s2 / s0 - (s1 / s0) ** 2).sum(axis=0)
        step = np.divide(score, information, out=np.zeros_like(score), where=information > 1e-12)
        step = np.clip(step, -2.0, 2.0)
        beta += step
        if np.nanmax(np.abs(step[valid])) < tolerance:
            break
    eta = np.clip(x * beta, -50.0, 50.0)
    s0 = np.cumsum(np.exp(eta), axis=0)[risk_rows]
    log_likelihood = eta[event_rows].sum(axis=0) - np.log(s0).sum(axis=0)
    standard_error = np.sqrt(np.divide(1.0, information, out=np.full_like(information, np.nan), where=information > 0))
    z = beta / standard_error
    p_value = 2.0 * norm.sf(np.abs(z))
    p_value[~valid] = np.nan
    log_likelihood[~valid] = np.nan
    beta[~valid] = np.nan
    return beta, standard_error, p_value, log_likelihood


def rank_block(
    x: pd.DataFrame,
    columns: list[str],
    y: pd.DataFrame,
    chunk_size: int,
    max_iter: int,
    tolerance: float,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    started = time.perf_counter()
    for start in range(0, len(columns), chunk_size):
        names = columns[start : start + chunk_size]
        values = x[names].to_numpy(dtype=np.float64, copy=False)
        beta, se, p_value, log_likelihood = univariate_cox_chunk(
            values,
            y["time"].to_numpy(dtype=np.float64),
            y["event"].to_numpy(dtype=np.int8),
            max_iter,
            tolerance,
        )
        rows.append(
            pd.DataFrame(
                {
                    "feature": names,
                    "coefficient": beta,
                    "standard_error": se,
                    "p_value": p_value,
                    "partial_log_likelihood": log_likelihood,
                }
            )
        )
        done = min(start + chunk_size, len(columns))
        print(f"  Cox taraması: {done}/{len(columns)} özellik", flush=True)
    result = pd.concat(rows, ignore_index=True)
    result = result.sort_values(
        ["p_value", "partial_log_likelihood"], ascending=[True, False], na_position="last"
    ).reset_index(drop=True)
    result["rank"] = np.arange(1, len(result) + 1)
    print(f"  Blok süresi: {time.perf_counter() - started:.2f} sn")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--max-iter", type=int, default=50)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.input_dir.resolve()
    started = time.perf_counter()
    x_train = read_frame(root / "X_train.parquet")
    x_test = read_frame(root / "X_test.parquet")
    x_train = collapse_ajcc_onehot(x_train, print_categories=True)
    x_test = collapse_ajcc_onehot(x_test, print_categories=False)
    y_train = read_frame(root / "y_train.parquet").loc[x_train.index]
    print(f"X_train başlangıç shape: {x_train.shape}; X_test: {x_test.shape}")
    print(f"Train örnek: {len(y_train)}; event: {int(y_train['event'].sum())}")
    clinical = [c for c in x_train.columns if c.startswith("clin_")]
    selected_by_block: dict[str, list[str]] = {}
    stats: dict[str, object] = {}
    for block, prefix in BLOCKS.items():
        columns = [c for c in x_train.columns if c.startswith(prefix)]
        print(f"{block}: {len(columns)} özellik için univariate CoxPH başlıyor")
        ranking = rank_block(
            x_train, columns, y_train, args.chunk_size, args.max_iter, args.tolerance
        )
        ranking.to_csv(root / f"univariate_cox_{block}.tsv", sep="\t", index=False)
        valid = ranking.dropna(subset=["p_value"])
        chosen_n = min(args.top_n, len(valid))
        chosen = valid.head(chosen_n)["feature"].tolist()
        selected_by_block[block] = chosen
        cutoff = float(valid.iloc[chosen_n - 1]["p_value"]) if chosen_n else None
        stats[block] = {
            "candidate_count": len(columns),
            "successfully_fitted": len(valid),
            "failed_or_constant": len(columns) - len(valid),
            "selected_count": chosen_n,
            "selection_rule": f"top-{args.top_n} by ascending Wald p-value; partial log-likelihood tie-break",
            "selected_boundary_p_value": cutoff,
        }
        print(
            f"{block}: {len(columns)} adaydan {chosen_n} seçildi; "
            f"top-N sınır p-value={cutoff}; başarısız/sabit={len(columns)-len(valid)}"
        )
    final_columns = [
        *selected_by_block["rna"],
        *selected_by_block["mirna"],
        *selected_by_block["methylation"],
        *clinical,
    ]
    selected_train = x_train[final_columns]
    selected_test = x_test[final_columns]
    write_frame(selected_train, root / "X_train_selected.parquet")
    write_frame(selected_test, root / "X_test_selected.parquet")
    payload = {
        "selection_scope": "train_only",
        "ranking": "univariate CoxPH Wald p-value ascending; partial log-likelihood tie-break",
        "top_n_per_omics": args.top_n,
        "blocks": stats,
        "selected_by_block": selected_by_block,
        "clinical_features_retained_without_univariate_filter": clinical,
        "final_features": final_columns,
        "final_feature_count": len(final_columns),
        "train_shape": list(selected_train.shape),
        "test_shape": list(selected_test.shape),
        "elapsed_seconds": time.perf_counter() - started,
    }
    (root / "selected_features.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Final X_train_selected shape: {selected_train.shape}")
    print(f"Final X_test_selected shape: {selected_test.shape}")
    print(f"Toplam seçim süresi: {payload['elapsed_seconds']:.2f} sn")


if __name__ == "__main__":
    main()
