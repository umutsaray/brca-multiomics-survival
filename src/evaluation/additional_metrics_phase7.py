#!/usr/bin/env python3
"""Recompute Table 9 metrics after replacing only DeepSurv with Phase 7 OOF risk."""
from __future__ import annotations
import json, sys, math
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from additional_metrics_v2 import evaluate, H, GRID  # noqa: E402
from foldstratified_test import corrected_oof  # noqa: E402


def main() -> None:
    repo = HERE.parents[1]
    outdir = repo / "results" / "phase7"; outdir.mkdir(parents=True, exist_ok=True)
    new = pd.read_parquet(outdir / "phase7_deepsurv_oof.parquet")
    old_paths = {
        "mixed": repo / "results/mixed_cohort/ablation_oof_extended.parquet",
        "primary": repo / "results/primary_tumor_only/ablation_oof_primary.parquet",
    }
    result = {"method": {
        "deep_surv_source": "Phase 7 leakage-corrected nested OOF predictions",
        "coxph_rsf_source": "retained Phase 2/5 OOF predictions; values must remain unchanged",
        "time_dependent_auc_horizons_days": H.tolist(),
        "ibs_grid": f"{len(GRID)} equally spaced points from day {GRID[0]:g} through day {GRID[-1]:g}",
        "ipcw_censoring_source": "other four OOF folds for each held-out fold",
    }}
    corrected_frames = {}
    for cohort, path in old_paths.items():
        data = corrected_oof(pd.read_parquet(path), new, cohort)
        corrected_frames[cohort] = data
        result[cohort] = {}
        for config in ("Clinical-only", "Full"):
            frame = data[data.configuration.eq(config)].reset_index(drop=True)
            print(f"{cohort} {config}: {len(frame)} rows")
            result[cohort][config] = evaluate(frame)

    prior = json.loads((repo / "results/primary_tumor_only/additional_metrics_v2.json").read_text(encoding="utf-8"))
    checks = []
    keys = ("mean_uno_c_index", "mean_auc_1y", "mean_auc_3y", "mean_auc_5y", "mean_weighted_mean_auc", "mean_ibs")
    for cohort in old_paths:
        for config in ("Clinical-only", "Full"):
            for model in ("coxph", "rsf"):
                for key in keys:
                    a = result[cohort][config][model][key]; b = prior[cohort][config][model][key]
                    if math.isnan(a) and math.isnan(b): diff = 0.0
                    elif math.isfinite(a) and math.isfinite(b): diff = abs(a-b)
                    else: diff = float("inf")
                    checks.append({"cohort": cohort, "configuration": config, "model": model, "metric": key,
                                   "phase7": a, "prior": b, "absolute_difference": diff})
    maxdiff = max(x["absolute_difference"] for x in checks)
    result["coxph_rsf_unchanged_verification"] = {"maximum_absolute_difference": maxdiff, "checks": checks}
    if maxdiff > 1e-12:
        raise RuntimeError(f"CoxPH/RSF changed unexpectedly; max absolute difference={maxdiff}")

    (outdir / "additional_metrics_phase7.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    lines = ["# Phase 7 updated additional metrics", "",
             "Only DeepSurv OOF risks were replaced. CoxPH and RSF were retained and verified exactly against the preceding results.", "",
             f"Maximum absolute CoxPH/RSF metric difference: `{maxdiff}`.", ""]
    for cohort in old_paths:
        lines += [f"## {cohort}", "", "| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |", "|---|---|---:|---:|---:|---:|---:|---:|"]
        for config in ("Clinical-only", "Full"):
            for model in ("coxph", "rsf", "deepsurv"):
                x = result[cohort][config][model]
                lines.append(f"| {config} | {model} | {x['mean_uno_c_index']:.4f} | {x['mean_auc_1y']:.4f} | {x['mean_auc_3y']:.4f} | {x['mean_auc_5y']:.4f} | {x['mean_weighted_mean_auc']:.4f} | {x['mean_ibs']:.4f} |")
        lines += [""]
    (outdir / "additional_metrics_phase7.md").write_text("\n".join(lines), encoding="utf-8")
    print("Wrote Phase 7 additional metrics; CoxPH/RSF unchanged verification passed")


if __name__ == "__main__": main()
