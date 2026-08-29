#!/usr/bin/env python3
"""Participant-level, fold-stratified paired bootstrap for Clinical-only vs Full."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored

MODELS = ("coxph", "rsf", "deepsurv")
CONFIGS = ("Clinical-only", "Full")


def cindex(frame: pd.DataFrame, model: str) -> float:
    return cindex_arrays(frame["event"].to_numpy(bool), frame["time"].to_numpy(float), frame[model].to_numpy(float))


def cindex_arrays(event: np.ndarray, time: np.ndarray, risk: np.ndarray) -> float:
    """Harrell C using the same implementation as all retained model outputs."""
    return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])


def corrected_oof(old: pd.DataFrame, new_dl: pd.DataFrame, cohort: str) -> pd.DataFrame:
    out = old[old["configuration"].isin(CONFIGS)].copy()
    repl = new_dl[new_dl["cohort"].eq(cohort)][
        ["sample_key", "fold", "configuration", "deepsurv"]
    ].rename(columns={"deepsurv": "deepsurv_phase7"})
    out = out.merge(repl, on=["sample_key", "fold", "configuration"], how="left", validate="one_to_one")
    if out["deepsurv_phase7"].isna().any():
        raise RuntimeError(f"Missing corrected DeepSurv predictions for {cohort}")
    out["deepsurv_old"] = out["deepsurv"]
    out["deepsurv"] = out.pop("deepsurv_phase7")
    out["participant"] = out["sample_key"].str[:12]
    return out


def fold_values(data: pd.DataFrame, model: str) -> dict:
    values = {}
    for fold in sorted(data["fold"].unique()):
        fd = data[data["fold"].eq(fold)]
        cc = cindex(fd[fd["configuration"].eq("Clinical-only")], model)
        cf = cindex(fd[fd["configuration"].eq("Full")], model)
        values[str(int(fold))] = {"clinical": cc, "full": cf, "delta_clinical_minus_full": cc - cf}
    deltas = np.array([v["delta_clinical_minus_full"] for v in values.values()])
    return {
        "folds": values,
        "mean_delta": float(deltas.mean()),
        "sample_sd_delta": float(deltas.std(ddof=1)),
        "mean_clinical": float(np.mean([v["clinical"] for v in values.values()])),
        "mean_full": float(np.mean([v["full"] for v in values.values()])),
    }


def bootstrap(data: pd.DataFrame, iterations: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    arrays = {m: np.empty(iterations, dtype=float) for m in MODELS}
    folds = sorted(data["fold"].unique())
    prepared = {}
    for fold in folds:
        fd = data[data["fold"].eq(fold)]
        by_config = {}
        for config in CONFIGS:
            cf = fd[fd["configuration"].eq(config)].sort_values("sample_key").reset_index(drop=True)
            groups = cf.groupby("participant", sort=True).indices
            participants = np.array(list(groups))
            by_config[config] = {
                "participants": participants,
                "groups": [np.asarray(groups[p], dtype=int) for p in participants],
                "event": cf["event"].to_numpy(bool), "time": cf["time"].to_numpy(float),
                **{m: cf[m].to_numpy(float) for m in MODELS},
            }
        if not np.array_equal(by_config[CONFIGS[0]]["participants"], by_config[CONFIGS[1]]["participants"]):
            raise RuntimeError(f"Configuration participant mismatch in fold {fold}")
        prepared[fold] = by_config
    for iteration in range(iterations):
        fold_deltas = {m: [] for m in MODELS}
        for fold in folds:
            by_config = prepared[fold]
            n_participants = len(by_config[CONFIGS[0]]["participants"])
            draw = rng.integers(0, n_participants, size=n_participants)
            for model in MODELS:
                scores = {}
                for config in CONFIGS:
                    item = by_config[config]
                    idx = np.concatenate([item["groups"][i] for i in draw])
                    scores[config] = cindex_arrays(item["event"][idx], item["time"][idx], item[model][idx])
                fold_deltas[model].append(scores["Clinical-only"] - scores["Full"])
        for model in MODELS:
            arrays[model][iteration] = np.mean(fold_deltas[model])
        if (iteration + 1) % 1000 == 0:
            print(f"  bootstrap {iteration + 1}/{iterations}", flush=True)
    result = {}
    for model, values in arrays.items():
        lo, hi = np.percentile(values, [2.5, 97.5])
        p_lo = (np.count_nonzero(values <= 0) + 1) / (iterations + 1)
        p_hi = (np.count_nonzero(values >= 0) + 1) / (iterations + 1)
        result[model] = {
            "ci_95": [float(lo), float(hi)],
            "two_sided_p": float(min(1.0, 2 * min(p_lo, p_hi))),
            "bootstrap_deltas": values.tolist(),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("results/phase7"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    new_dl = pd.read_parquet("results/phase7/phase7_deepsurv_oof.parquet")
    paths = {
        "mixed": "results/mixed_cohort/ablation_oof_extended.parquet",
        "primary": "results/primary_tumor_only/ablation_oof_primary.parquet",
    }
    payload = {
        "method": "participant bootstrap within each outer fold; mean of five fold-specific C-index deltas",
        "delta_direction": "Clinical-only minus Full",
        "iterations": args.iterations,
        "seed": args.seed,
        "cohorts": {},
    }
    md = [
        "# Fold-stratified clinical-only versus full multi-omics comparison",
        "",
        f"Paired participant-level bootstrap within each outer fold ({args.iterations:,} iterations; seed {args.seed}). "
        "The estimand is the unweighted mean of the five fold-specific C-index differences: Clinical-only minus Full. "
        "Positive values favor Clinical-only.", "",
    ]
    for cohort, path in paths.items():
        print(f"[{cohort}] loading OOF predictions", flush=True)
        old = pd.read_parquet(path)
        data = corrected_oof(old, new_dl, cohort)
        summaries = {model: fold_values(data, model) for model in MODELS}
        print(f"[{cohort}] running {args.iterations} paired fold-stratified bootstrap iterations", flush=True)
        boot = bootstrap(data, args.iterations, args.seed)
        for model in MODELS:
            summaries[model].update({k: v for k, v in boot[model].items() if k != "bootstrap_deltas"})
        payload["cohorts"][cohort] = {
            "n_rows": int(data[data["configuration"].eq("Full")].shape[0]),
            "n_participants": int(data["participant"].nunique()),
            "models": summaries,
            "bootstrap_deltas": {m: boot[m]["bootstrap_deltas"] for m in MODELS},
        }
        md += [f"## {cohort}", "", "| Model | Mean Clinical C | Mean Full C | Mean delta | Fold SD | 95% CI | Two-sided p | Interpretation |", "|---|---:|---:|---:|---:|---:|---:|---|"]
        for model in MODELS:
            s = summaries[model]; lo, hi = s["ci_95"]
            interpretation = "Clinical-only favored" if lo > 0 else ("Full favored" if hi < 0 else "No statistically detectable difference")
            md.append(f"| {model} | {s['mean_clinical']:.4f} | {s['mean_full']:.4f} | {s['mean_delta']:.4f} | {s['sample_sd_delta']:.4f} | [{lo:.4f}, {hi:.4f}] | {s['two_sided_p']:.4g} | {interpretation} |")
        md += ["", "Fold-level values:", "", "| Model | Fold | Clinical C | Full C | Delta |", "|---|---:|---:|---:|---:|"]
        for model in MODELS:
            for fold, row in summaries[model]["folds"].items():
                md.append(f"| {model} | {fold} | {row['clinical']:.4f} | {row['full']:.4f} | {row['delta_clinical_minus_full']:.4f} |")
        md.append("")

    json_path = args.output_dir / "foldstratified_test_results.json"
    md_path = args.output_dir / "foldstratified_test_results.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md_path.write_text("\n".join(md), encoding="utf-8")
    print(f"Wrote {json_path} and {md_path}")


if __name__ == "__main__":
    main()
