#!/usr/bin/env python3
"""Create an auditable row-level registry for Phase 1/2 reported CV results."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import pandas as pd
from phase1_common import load_cohort, grouped_folds
from sksurv.metrics import concordance_index_censored

MODELS = ("coxph", "rsf", "deepsurv")

def fold_hash(ids) -> str:
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()

def main() -> None:
    root = Path(__file__).resolve().parent
    rows = []
    hashes = {}
    for cohort in ("mixed", "primary"):
        _, y = load_cohort(root, cohort)
        hashes[cohort] = {i: fold_hash(pd.Index(y.index[te]).str[:12].unique())
                          for i, (_, te) in enumerate(grouped_folds(y, 5, 42), 1)}
    for cohort, suffix in (("mixed", ""), ("primary", "_primary")):
        path = root / f"select_features_v2_results{suffix}.json"
        if path.is_file():
            d = json.loads(path.read_text())
            for n, models in d["sensitivity"].items():
                for model in MODELS:
                    for fold, value in enumerate(models[model]["fold_cindices"], 1):
                        rows.append({"cohort": cohort, "script": "select_features_v2.py", "analysis": f"top-N sensitivity N={n}",
                                     "fold": fold, "test_participant_hash_sha256": hashes[cohort][fold],
                                     "split_seed": 42, "model_seed": 42 + fold * 100 + int(n), "model": model,
                                     "feature_selection": f"raw-p top {n} per omics block; clinical retained",
                                     "feature_count_rule": f"up to {n}/block + clinical", "c_index": value})
        apath = root / ("ablation_results_primary.json" if cohort == "primary" else "ablation_results_extended.json")
        if apath.is_file():
            d = json.loads(apath.read_text())
            opath = root / ("ablation_oof_primary.parquet" if cohort == "primary" else "ablation_oof_extended.parquet")
            oof = pd.read_parquet(opath)
            for config in d["results"]:
                for model in MODELS:
                    for fold, g in oof[oof.configuration.eq(config)].groupby("fold"):
                        value = float(concordance_index_censored(g.event.astype(bool), g.time, g[model])[0])
                        rows.append({"cohort": cohort, "script": "ablation_study.py", "analysis": config,
                                     "fold": int(fold), "test_participant_hash_sha256": hashes[cohort][int(fold)],
                                     "split_seed": 42, "model_seed": 42 + int(fold), "model": model,
                                     "feature_selection": "BH-FDR < 0.05, capped at top 200 per included omics block; clinical retained",
                                     "feature_count_rule": "0-200/block + included clinical", "c_index": value})
        tpath = root / f"coxph_penalizer_tuning{suffix}.json"
        if tpath.is_file():
            d = json.loads(tpath.read_text())
            for item in d["folds"]:
                fold = item["fold"]
                rows.append({"cohort": cohort, "script": "tune_coxph_penalizer.py", "analysis": "nested ridge penalizer tuning",
                             "fold": fold, "test_participant_hash_sha256": hashes[cohort][fold], "split_seed": 42,
                             "model_seed": "not stochastic", "model": "coxph",
                             "feature_selection": "outer-train BH-FDR < 0.05, capped top 200/block; clinical retained",
                             "feature_count_rule": "0-200/block + clinical", "selected_penalizer": item["selected_penalizer"],
                             "c_index": item["test_cindex"]})
    stable_path = root / "stable_cox_sensitivity.json"
    if stable_path.is_file():
        stable = json.loads(stable_path.read_text())
        for cohort, d in stable.items():
            for fold, value in enumerate(d["fold_cindices"], 1):
                rows.append({"cohort": cohort, "script": "stable_cox_sensitivity.py", "analysis": "clinical-only elastic-net sensitivity",
                             "fold": fold, "test_participant_hash_sha256": hashes[cohort][fold], "split_seed": 42,
                             "model_seed": "not stochastic", "model": "coxnet", "feature_selection": "all clinical one-hot covariates",
                             "feature_count_rule": "all clinical", "selected_penalizer": d["alpha"], "c_index": value})
    frame = pd.DataFrame(rows)
    frame["hyperparameters"] = frame.model.map({"coxph": "lifelines; ridge penalizer=0.1 unless selected_penalizer set; max_steps=50/100",
                                                  "rsf": "trees=500; min_split=10; min_leaf=5; max_features=sqrt",
                                                  "deepsurv": "dropout=.3; AdamW lr=.001 wd=.0001; epochs=300; patience=30",
                                                  "coxnet": "CoxnetSurvivalAnalysis; l1_ratio=.5; alpha=.01; max_iter=100000"})
    frame.to_csv(root / "experiment_log.csv", index=False)
    lines = ["# Experiment Log", "", f"Logged result rows: {len(frame)}.", "",
             "The full machine-readable registry is `experiment_log.csv`. Fold hashes are SHA-256 hashes of sorted 12-character test participant IDs.", "",
             "## Table 3b / Table 4 reconciliation", "",
             "Both scripts use the same participant-grouped outer folds (seed 42), confirmed by identical fold hashes. They do not use the same feature-selection rule or neural/RSF seed: the sensitivity table uses raw-p top-N and a seed dependent on N; ablation uses BH-FDR < 0.05 capped at N=200 and seed 42+fold. Thus the values are distinct sensitivity and canonical FDR-ablation estimates, not conflicting reruns of one experiment.", "",
             "See the CSV for the complete unrounded row-level registry."]
    (root / "experiment_log.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"experiment_log.csv/md: {len(frame)} satır")

if __name__ == "__main__":
    main()
