#!/usr/bin/env python3
"""Compile Phase 7 DeepSurv, bootstrap, SHAP, and provenance reports."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
from foldstratified_test import corrected_oof, cindex  # noqa: E402


def f4(x): return f"{x:.4f}"


def main() -> None:
    out = REPO / "results/phase7"; out.mkdir(parents=True, exist_ok=True)
    new_dl = pd.read_parquet(out / "phase7_deepsurv_oof.parquet")
    phase7 = json.loads((out / "phase7_deepsurv_results.json").read_text(encoding="utf-8"))
    strat = json.loads((out / "foldstratified_test_results.json").read_text(encoding="utf-8"))
    pooled = json.loads((REPO / "results/primary_tumor_only/bootstrap_9999_primary.json").read_text(encoding="utf-8"))
    old_shap = json.loads((REPO / "results/primary_tumor_only/shap_v2_results.json").read_text(encoding="utf-8"))
    new_shap = json.loads((out / "shap_phase7_results.json").read_text(encoding="utf-8"))
    old_paths = {
        "mixed": REPO / "results/mixed_cohort/ablation_oof_extended.parquet",
        "primary": REPO / "results/primary_tumor_only/ablation_oof_primary.parquet",
    }

    # Expand fold-stratified report with the required old-pooled versus new table.
    fs = ["# Fold-stratified clinical-only versus full multi-omics comparison", "",
          "The new analysis resamples participants independently within every outer fold, never compares raw risks across folds, computes Harrell C within each fold, and uses the unweighted mean of five fold-specific deltas. Delta is Clinical-only minus Full; positive values favor Clinical-only. All models share the same participant resample indices in each of 9,999 iterations (seed 42).", "",
          "## Old pooled versus new fold-stratified results", "",
          "The old stored estimates used Full minus Clinical. They are sign-normalized below to Clinical minus Full for a direct comparison.", "",
          "| Cohort | Model | Old pooled delta | Old 95% CI | Old p | New fold-stratified delta | New 95% CI | New p |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for cohort in ("mixed", "primary"):
        old_by_model = {x["model"]: x for x in pooled[cohort]["comparisons"]}
        for model in ("coxph", "rsf", "deepsurv"):
            o = old_by_model[model]; n = strat["cohorts"][cohort]["models"][model]
            old_delta = -o["delta_c_index"]; old_ci = [-o["ci_95"][1], -o["ci_95"][0]]
            fs.append(f"| {cohort} | {model} | {f4(old_delta)} | [{f4(old_ci[0])}, {f4(old_ci[1])}] | {o['bootstrap_p_value']:.4g} | {f4(n['mean_delta'])} | [{f4(n['ci_95'][0])}, {f4(n['ci_95'][1])}] | {n['two_sided_p']:.4g} |")
    fs += ["", "## Fold-level new results", ""]
    for cohort in ("mixed", "primary"):
        fs += [f"### {cohort}", "", "| Model | Fold | Clinical C | Full C | Delta |", "|---|---:|---:|---:|---:|"]
        for model in ("coxph", "rsf", "deepsurv"):
            s = strat["cohorts"][cohort]["models"][model]
            for fold, row in s["folds"].items():
                fs.append(f"| {model} | {fold} | {f4(row['clinical'])} | {f4(row['full'])} | {f4(row['delta_clinical_minus_full'])} |")
        fs.append("")
    fs += ["## Direction and significance", ""]
    for cohort in ("mixed", "primary"):
        for model in ("coxph", "rsf", "deepsurv"):
            n = strat["cohorts"][cohort]["models"][model]; lo, hi = n["ci_95"]
            verdict = "Clinical-only is favored and the CI excludes zero" if lo > 0 else ("Full is favored and the CI excludes zero" if hi < 0 else "the CI includes zero; no statistically detectable difference")
            fs.append(f"- {cohort}, {model}: delta {f4(n['mean_delta'])}, 95% CI [{f4(lo)}, {f4(hi)}], p={n['two_sided_p']:.4g}; {verdict}.")
    fs += ["", "The prior pooled values are retained only for audit comparison; they are not the Phase 7 primary estimates."]
    (out / "foldstratified_test_results.md").write_text("\n".join(fs), encoding="utf-8")

    # Leakage-fix report: every old/new DeepSurv fold and the SHAP comparison.
    lines = ["# DeepSurv internal-validation leakage correction", "",
             "Phase 7 splits each outer-training set into participant-grouped, event-stratified internal-train and internal-validation subsets before any learned preprocessing or outcome-guided feature selection. Internal validation is used only for early stopping; each outer test fold is evaluated once.", "",
             "Scope: exactly 20 retrainings (2 cohorts × Clinical-only/Full × 5 folds). Other ablation configurations were not rerun.", "",
             "## Old versus new DeepSurv C-index", "", "| Cohort | Configuration | Fold | Old C | New C | New − old |", "|---|---|---:|---:|---:|---:|"]
    summary_rows = []
    for cohort, path in old_paths.items():
        old = pd.read_parquet(path)
        corrected = corrected_oof(old, new_dl, cohort)
        for config in ("Clinical-only", "Full"):
            part = corrected[corrected.configuration.eq(config)].copy()
            old_cs=[]; new_cs=[]
            for fold in map(str, range(1,6)):
                fd=part[part.fold.eq(int(fold))]
                old_fd=fd.copy();old_fd["deepsurv"]=old_fd["deepsurv_old"]
                oc=cindex(old_fd,"deepsurv");nc=cindex(fd,"deepsurv")
                old_cs.append(oc);new_cs.append(nc)
                lines.append(f"| {cohort} | {config} | {fold} | {f4(oc)} | {f4(nc)} | {f4(nc-oc)} |")
            summary_rows.append((cohort,config,np.mean(old_cs),np.std(old_cs,ddof=1),np.mean(new_cs),np.std(new_cs,ddof=1)))
    lines += ["", "Mean ± sample SD:", "", "| Cohort | Configuration | Old | New |", "|---|---|---:|---:|"]
    for cohort,config,om,os,nm,ns in summary_rows:
        lines.append(f"| {cohort} | {config} | {f4(om)} ± {f4(os)} | {f4(nm)} ± {f4(ns)} |")

    old_features=[x["feature"] for x in old_shap["aggregate_top20"]];new_features=[x["feature"] for x in new_shap["aggregate_top20"]]
    overlap=len(set(old_features)&set(new_features)); union=len(set(old_features)|set(new_features)); old_new_j=overlap/union
    lines += ["", "## SHAP old versus new", "",
              f"Old cross-fold mean pairwise top-20 Jaccard: `{old_shap['mean_pairwise_jaccard']}`.",
              f"New cross-fold mean pairwise top-20 Jaccard: `{new_shap['mean_pairwise_jaccard']}`.",
              f"Old/new aggregate top-20 overlap: `{overlap}/20`; Jaccard: `{old_new_j}`.", "",
              "| Rank | Old aggregate feature | New aggregate feature | New mean |SHAP| | New folds selected |", "|---:|---|---|---:|---:|"]
    for i,(o,n) in enumerate(zip(old_shap["aggregate_top20"],new_shap["aggregate_top20"]),1):
        lines.append(f"| {i} | {o['feature']} | {n['feature']} | {n['mean_abs_shap_across_folds_where_selected']:.8g} | {n['folds_selected']}/5 |")

    lines += ["", "## Primary finding after correction", ""]
    for cohort in ("mixed","primary"):
        d=strat["cohorts"][cohort]["models"]["deepsurv"]; lo,hi=d["ci_95"]
        conclusion="Clinical-only superiority retained" if lo>0 else ("Full superiority" if hi<0 else "direction may favor Clinical-only, but statistical significance is not retained")
        lines.append(f"- {cohort}: Clinical − Full = {f4(d['mean_delta'])} (95% CI [{f4(lo)}, {f4(hi)}], p={d['two_sided_p']:.4g}); {conclusion}.")
    lines += ["", "## Seeds and split schema", "",
              "- Outer StratifiedGroupKFold: seed 42.",
              "- Internal 7-fold StratifiedGroupKFold first held split: seed = 4200 + outer-fold (4201–4205).",
              "- DeepSurv initialization/training: seed = 42 + outer-fold (43–47).",
              "- Fold-stratified paired participant bootstrap: seed 42, 9,999 iterations.",
              "", "CoxPH and RSF were not retrained. Their Table 9 metrics were recomputed from retained OOF predictions and matched the prior results with maximum absolute difference 0.0."]
    (out / "deepsurv_leakage_fix_results.md").write_text("\n".join(lines), encoding="utf-8")

    # Current SHAP table and supplementary current ranking.
    shap_md=["# Phase 7 leakage-corrected SHAP aggregate top 20", "", f"Mean pairwise fold top-20 Jaccard: `{new_shap['mean_pairwise_jaccard']}`.", "", "| Rank | Feature | Mean |SHAP| across folds where selected | Folds selected |", "|---:|---|---:|---:|"]
    csv=[]
    for rank,row in enumerate(new_shap["aggregate_top20"],1):
        shap_md.append(f"| {rank} | {row['feature']} | {row['mean_abs_shap_across_folds_where_selected']:.8g} | {row['folds_selected']}/5 |")
        csv.append({"rank":rank,**row})
    (out / "shap_phase7_results.md").write_text("\n".join(shap_md),encoding="utf-8")
    pd.DataFrame(csv).to_csv(REPO / "supplementary/table_s1_shap_full_ranking.csv",index=False)
    print("Compiled Phase 7 reports and current SHAP table")


if __name__ == "__main__": main()
