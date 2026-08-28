#!/usr/bin/env python3
"""Diagnose fold-wise OOF risk scale offsets without inventing unavailable data."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored

MODELS = ("coxph", "rsf", "deepsurv")
DEFINITION = {"coxph": "lifelines predict_partial_hazard (partial hazard)",
              "rsf": "scikit-survival RandomSurvivalForest.predict risk score",
              "deepsurv": "raw neural-network scalar output (log-risk)"}

def ci(d, risk):
    return float(concordance_index_censored(d.event.astype(bool), d.time, risk)[0])

def analyze(path: Path) -> dict:
    d = pd.read_parquet(path)
    result = {"source": path.name, "n_rows": len(d), "scores": {}}
    for model in MODELS:
        fold_stats = []
        centered = d[model].copy()
        for fold, g in d.groupby("fold"):
            fold_stats.append({"fold": int(fold), "n": len(g), "mean": float(g[model].mean()),
                               "sample_variance": float(g[model].var(ddof=1))})
            centered.loc[g.index] = g[model] - g[model].mean()
        result["scores"][model] = {"definition": DEFINITION[model], "fold_statistics": fold_stats,
                                    "pooled_c_index_uncentered": ci(d, d[model].to_numpy()),
                                    "pooled_c_index_test_fold_mean_centered_diagnostic_only": ci(d, centered.to_numpy())}
    return result

def main():
    root = Path(__file__).resolve().parent
    payload = {}
    for cohort, suffix in (("mixed", ""), ("primary", "_primary")):
        path = root / f"phase1_oof_predictions{suffix}.parquet"
        if path.is_file(): payload[cohort] = analyze(path)
    payload["methodological_note"] = ("The completed OOF files do not store training-fold risk means. "
        "Therefore a leakage-safe training-mean-centered pooled C-index cannot be reconstructed and is not fabricated here. "
        "The supplied centered value subtracts each held-out fold's own mean and is labelled diagnostic-only; it must not be used as a performance estimate. "
        "Future model fits must persist the training-risk mean and subtract that constant from their held-out risks.")
    (root / "oof_scale_diagnostics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = ["# OOF Risk-scale Diagnostics", "", payload["methodological_note"], ""]
    for cohort, d in payload.items():
        if cohort == "methodological_note": continue
        lines += [f"## {cohort}", "", "| Model | Fold | N | Mean | Sample variance |", "|---|---:|---:|---:|---:|"]
        for model, x in d["scores"].items():
            for f in x["fold_statistics"]: lines.append(f"| {model} | {f['fold']} | {f['n']} | {f['mean']} | {f['sample_variance']} |")
        lines += ["", "| Model | Risk definition | Pooled uncentered C | Test-fold-mean centered C (diagnostic only) |", "|---|---|---:|---:|"]
        for model, x in d["scores"].items(): lines.append(f"| {model} | {x['definition']} | {x['pooled_c_index_uncentered']} | {x['pooled_c_index_test_fold_mean_centered_diagnostic_only']} |")
        lines.append("")
    (root / "oof_scale_diagnostics.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(payload, indent=2))

if __name__ == "__main__": main()
