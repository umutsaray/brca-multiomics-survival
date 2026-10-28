#!/usr/bin/env python3
"""Phase 17 Task 6: reduced CoxPH-only training-size learning curve."""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import lifelines
import numpy as np
import pandas as pd
import sklearn
import sksurv
from lifelines import CoxPHFitter
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sksurv.metrics import concordance_index_censored

HERE = Path(__file__).resolve().parent
for sub in ("common", "data_prep", "preprocessing", "feature_selection", "models", "evaluation"):
    sys.path.insert(0, str(HERE.parents[0] / sub))
from phase1_common import grouped_folds, load_cohort, rank_fold, select_from_rankings  # noqa: E402
from phase15_task2_receptor_sensitivity import prepare_outer_uncached  # noqa: E402

CONFIGS = ("Clinical-only", "Full")
FRACTIONS = ((1, 0.50), (2, 0.75))
REPETITIONS = (1, 2, 3)
BLOCK_PREFIXES = {"rna": "rna_", "mirna": "mirna_", "methylation": "meth_", "clinical": "clin_"}


def cindex(event, duration, risk):
    return float(concordance_index_censored(np.asarray(event, bool), np.asarray(duration, float), np.asarray(risk, float))[0])


def participant_subsample(y, outer_train, fraction, seed):
    frame = y.iloc[outer_train].copy()
    frame["participant"] = frame.index.astype(str).str[:12]
    participant = frame.groupby("participant", sort=True).agg(event=("event", "max"))
    chosen, _ = train_test_split(
        participant.index.to_numpy(), train_size=fraction, random_state=seed,
        stratify=participant.event.to_numpy(), shuffle=True,
    )
    chosen = set(chosen)
    index = np.asarray([i for i in outer_train if str(y.index[i])[:12] in chosen], dtype=int)
    return index, sorted(chosen)


def fit_cox(xfit, yfit, xtest, ytest, penalizer, max_steps):
    event = yfit.event.astype(bool)
    event_var = xfit.loc[event].var()
    censored_var = xfit.loc[~event].var()
    columns = [c for c in xfit if event_var[c] > 1e-12 and censored_var[c] > 1e-12]
    scaler = StandardScaler()
    train = scaler.fit_transform(xfit[columns])
    test = scaler.transform(xtest[columns])
    fit = pd.DataFrame(train, columns=columns)
    fit["time"] = yfit.time.to_numpy()
    fit["event"] = yfit.event.to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = CoxPHFitter(penalizer=penalizer)
        model.fit(fit, duration_col="time", event_col="event", fit_options={"max_steps": max_steps})
    risk = model.predict_partial_hazard(pd.DataFrame(test, columns=columns)).to_numpy().ravel()
    return {
        "risk": risk.tolist(), "cindex": cindex(ytest.event, ytest.time, risk),
        "columns_after_event_censor_variance_filter": len(columns),
        "warnings": [str(w.message) for w in caught],
    }


def counts(columns):
    return {name: int(sum(c.startswith(prefix) for c in columns)) for name, prefix in BLOCK_PREFIXES.items()}


def checkpoint_path(checkpoint_dir, fold, fraction, repetition):
    return checkpoint_dir / f"fold{fold}_fraction{int(fraction * 100)}_rep{repetition}.json"


def stored_100(repo, y, folds):
    oof = pd.read_parquet(repo / "results/phase9/ablation_oof_phase9.parquet")
    oof = oof[(oof.cohort == "primary") & oof.configuration.isin(CONFIGS)].copy()
    t3 = json.loads((repo / "results/phase16/task3_fusion_strategies.json").read_text(encoding="utf-8"))
    records = []
    for fold, (train_idx, test_idx) in enumerate(folds, 1):
        selected = t3["folds"][fold - 1]["baseline_fdr_selected_features"]
        for config in CONFIGS:
            config_selected = [c for c in selected if c.startswith("clin_")] if config == "Clinical-only" else selected
            selection_counts = counts(config_selected)
            d = oof[(oof.fold == fold) & (oof.configuration == config)].sort_values("sample_key")
            records.append({
                "fold": fold, "fraction": 1.0, "repetition": None, "seed": 42,
                "configuration": config, "training_participants": int(pd.Index(y.index[train_idx]).str[:12].nunique()),
                "training_rows": int(len(train_idx)), "training_events": int(y.iloc[train_idx].event.sum()),
                "selected_features_per_block": selection_counts,
                "selected_feature_total": len(config_selected), "test_sample_keys": d.sample_key.tolist(),
                "test_event": d.event.astype(int).tolist(), "test_time": d.time.astype(float).tolist(),
                "risk": d.coxph.astype(float).tolist(), "cindex": cindex(d.event, d.time, d.coxph),
                "source": "stored Table 7 CoxPH OOF predictions",
            })
    return records


def fold_summary(records):
    result = {}
    for fraction in (0.5, 0.75, 1.0):
        result[str(fraction)] = {}
        for config in CONFIGS:
            fold_values, fold_detail = [], []
            for fold in range(1, 6):
                r = [x for x in records if x["fraction"] == fraction and x["configuration"] == config and x["fold"] == fold]
                if not r:
                    continue
                vals = [x["cindex"] for x in r]
                fold_values.append(float(np.mean(vals)))
                fold_detail.append({
                    "fold": fold, "repetition_cindices": vals, "repetition_mean_cindex": float(np.mean(vals)),
                    "training_participants": [x["training_participants"] for x in r],
                    "training_events": [x["training_events"] for x in r],
                    "selected_features_per_block": [x["selected_features_per_block"] for x in r],
                })
            result[str(fraction)][config] = {
                "completed_folds": len(fold_values), "fold_repetition_mean_cindices": fold_values,
                "mean": float(np.mean(fold_values)) if fold_values else None,
                "sample_sd": float(np.std(fold_values, ddof=1)) if len(fold_values) > 1 else None,
                "folds": fold_detail,
            }
    return result


def paired_bootstrap(records, fraction, iterations, seed):
    by = {(r["fold"], r["configuration"], r["repetition"]): r for r in records if r["fraction"] == fraction}
    expected_reps = [None] if fraction == 1.0 else list(REPETITIONS)
    if any((fold, config, rep) not in by for fold in range(1, 6) for config in CONFIGS for rep in expected_reps):
        return None
    rng = np.random.default_rng(seed)
    deltas = np.empty(iterations)
    for iteration in range(iterations):
        fold_deltas = []
        for fold in range(1, 6):
            ref = by[(fold, CONFIGS[0], expected_reps[0])]
            participants = np.asarray([str(k)[:12] for k in ref["test_sample_keys"]])
            unique = pd.unique(participants)
            chosen = unique[rng.integers(0, len(unique), len(unique))]
            indices = np.concatenate([np.flatnonzero(participants == p) for p in chosen])
            scores = {}
            for config in CONFIGS:
                rep_scores = []
                for rep in expected_reps:
                    rec = by[(fold, config, rep)]
                    rep_scores.append(cindex(np.asarray(rec["test_event"])[indices], np.asarray(rec["test_time"])[indices], np.asarray(rec["risk"])[indices]))
                scores[config] = float(np.mean(rep_scores))
            fold_deltas.append(scores["Clinical-only"] - scores["Full"])
        deltas[iteration] = np.mean(fold_deltas)
    point_folds=[]
    for fold in range(1,6):
        scores={config:float(np.mean([by[(fold,config,rep)]["cindex"] for rep in expected_reps])) for config in CONFIGS}
        point_folds.append(scores["Clinical-only"]-scores["Full"])
    return {
        "estimand": "unweighted mean across folds of repetition-averaged within-fold Harrell-C differences; Clinical-only minus Full",
        "iterations": iterations, "seed": seed, "point_difference":float(np.mean(point_folds)),
        "ci_95": [float(x) for x in np.percentile(deltas, [2.5, 97.5])],
        "bootstrap_deltas": deltas.tolist(),
    }


def trend_label(values):
    if any(v is None for v in values):
        return None
    if values[0] == values[1] == values[2]:
        return "flat"
    if values[0] <= values[1] <= values[2]:
        return "increasing"
    if values[0] >= values[1] >= values[2]:
        return "decreasing"
    return "non-monotonic"


def write_outputs(out, records, args, started, status, stop_reason=None, projected=None):
    summary = fold_summary(records)
    boots = {str(f): paired_bootstrap(records, f, args.bootstrap_iterations, args.bootstrap_seed) for f in (0.5, 0.75, 1.0)}
    trends = {c: {"means_50_75_100": [summary[str(f)][c]["mean"] for f in (0.5, 0.75, 1.0)]} for c in CONFIGS}
    for c in CONFIGS:
        trends[c]["label"] = trend_label(trends[c]["means_50_75_100"])
    payload = {
        "task": "Phase 17 Task 6 reduced learning curve", "status": status, "stop_reason": stop_reason,
        "projected_remaining_seconds_at_stop": projected,
        "design_fixed_before_run": {
            "cohort": "corrected primary-tumor-only", "models_run": ["CoxPH"],
            "models_skipped": {"RSF": "pre-fixed runtime decision", "DeepSurv": "pre-fixed runtime decision"},
            "fractions_rerun": [0.5, 0.75], "fraction_1_source": "stored Table 7",
            "repetitions": 3, "seed_formula": "7000 + 100*outer_fold + 10*fraction_index + repetition",
            "outer_folds": "same participant-grouped event-stratified 5 folds as Table 7; seed 42",
            "subsampling": "participant-level and event-stratified; outer test folds unchanged",
            "preprocessing_selection": "fit only on each subsampled outer-training fit subset",
            "cox_penalizer": args.cox_penalizer, "cox_max_steps": args.cox_max_steps,
            "top_n_per_omics_block": args.top_n, "FDR_BH_lt": 0.05,
            "bootstrap_iterations": args.bootstrap_iterations, "bootstrap_seed": args.bootstrap_seed,
            "runtime_stop_rule_seconds": args.max_projected_remaining_seconds,
        },
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "sklearn": sklearn.__version__, "sksurv": sksurv.__version__, "lifelines": lifelines.__version__},
        "command": " ".join([sys.executable, *sys.argv]), "runtime_seconds": time.perf_counter() - started,
        "completed_checkpoint_count": sum(r["fraction"] < 1 for r in records) // len(CONFIGS), "summary": summary,
        "clinical_minus_full_bootstrap": boots, "trends": trends,
    }
    (out / "task6_learning_curve.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = [
        "# Phase 17 Task 6 — Reduced CoxPH-only learning curve", "", f"Status: **{status}**.", "",
        "RSF and DeepSurv were skipped as a pre-fixed runtime decision. Fractions 50% and 75% use three repetitions; 100% is the stored Table 7 CoxPH result.", "",
    ]
    if stop_reason:
        lines += [f"Stop reason: {stop_reason}", f"Projected remaining seconds at stop: `{projected}`.", ""]
    lines += ["## C-index summary", "", "| Fraction | Configuration | Completed folds | Mean | Sample SD | Fold values |", "|---:|---|---:|---:|---:|---|"]
    for f in (0.5, 0.75, 1.0):
        for c in CONFIGS:
            z = summary[str(f)][c]
            lines.append(f"| {f} | {c} | {z['completed_folds']} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_repetition_mean_cindices'])}` |")
    lines += ["", "## Clinical-only minus Full bootstrap", "", "| Fraction | Point difference | 95% CI |", "|---:|---:|---|"]
    for f in (0.5, 0.75, 1.0):
        z = boots[str(f)]
        lines.append(f"| {f} | {None if z is None else z['point_difference']} | `{json.dumps(None if z is None else z['ci_95'])}` |")
    lines += ["", "## Observed trends", ""]
    for c, z in trends.items():
        lines.append(f"- {c}: `{z['label']}`; means at 50%, 75%, 100% = `{json.dumps(z['means_50_75_100'])}`.")
    lines += ["", "Per-fold training participant/event counts, repetition-level C-indices, selected feature counts per block, raw risks, and all bootstrap draws are retained without rounding in the JSON/checkpoints.", "", f"Actual elapsed time (seconds): `{payload['runtime_seconds']}`."]
    (out / "task6_learning_curve.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return payload


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--repo-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--outer-seed", type=int, default=42)
    p.add_argument("--meth-prefilter", type=int, default=10000)
    p.add_argument("--methylation-chunksize", type=int, default=5000)
    p.add_argument("--knn-neighbors", type=int, default=5)
    p.add_argument("--top-n", type=int, default=200)
    p.add_argument("--cox-penalizer", type=float, default=0.1)
    p.add_argument("--cox-max-steps", type=int, default=50)
    p.add_argument("--bootstrap-iterations", type=int, default=9999)
    p.add_argument("--bootstrap-seed", type=int, default=42)
    p.add_argument("--max-projected-remaining-seconds", type=float, default=10800.0)
    a = p.parse_args()
    started = time.perf_counter()
    root, repo, out = a.input_dir.resolve(), a.repo_dir.resolve(), a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    checkpoints = out / "task6_checkpoints"
    checkpoints.mkdir(exist_ok=True)
    base, y = load_cohort(root, "primary")
    folds = grouped_folds(y, 5, a.outer_seed)
    records = stored_100(repo, y, folds)
    for path in sorted(checkpoints.glob("fold*_fraction*_rep*.json")):
        records.extend(json.loads(path.read_text(encoding="utf-8"))["configurations"])
    completed = sum(r["fraction"] < 1 for r in records) // 2
    total_jobs = 5 * len(FRACTIONS) * len(REPETITIONS)
    prep_args = SimpleNamespace(meth_prefilter=a.meth_prefilter, methylation_chunksize=a.methylation_chunksize, knn_neighbors=a.knn_neighbors)
    for fold, (outer_train, outer_test) in enumerate(folds, 1):
        for fraction_index, fraction in FRACTIONS:
            for repetition in REPETITIONS:
                cp = checkpoint_path(checkpoints, fold, fraction, repetition)
                if cp.is_file():
                    continue
                job_started = time.perf_counter()
                seed = 7000 + 100 * fold + 10 * fraction_index + repetition
                fit_idx, participants = participant_subsample(y, outer_train, fraction, seed)
                xfit, xtest, preprocessing = prepare_outer_uncached(root, base, y, fit_idx, outer_test, prep_args)
                yfit, ytest = y.iloc[fit_idx], y.iloc[outer_test]
                rankings = rank_fold(xfit, yfit)
                selected, selection = select_from_rankings(rankings, xfit, a.top_n, True)
                configs = []
                for config in CONFIGS:
                    columns = [c for c in selected if c.startswith("clin_")] if config == "Clinical-only" else selected
                    result = fit_cox(xfit[columns], yfit, xtest[columns], ytest, a.cox_penalizer, a.cox_max_steps)
                    configs.append({
                        "fold": fold, "fraction": fraction, "fraction_index": fraction_index, "repetition": repetition,
                        "seed": seed, "configuration": config, "training_participants": len(participants),
                        "training_rows": int(len(fit_idx)), "training_events": int(yfit.event.sum()),
                        "selected_features_per_block": counts(columns), "selected_feature_total": len(columns),
                        "selection_diagnostics": selection, "preprocessing": preprocessing,
                        "test_sample_keys": ytest.index.tolist(), "test_event": ytest.event.astype(int).tolist(),
                        "test_time": ytest.time.astype(float).tolist(), **result,
                        "source": "Phase 17 fold-local subsample fit",
                    })
                checkpoint = {
                    "status": "completed", "fold": fold, "fraction": fraction, "repetition": repetition, "seed": seed,
                    "job_runtime_seconds": time.perf_counter() - job_started, "configurations": configs,
                }
                cp.write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")
                records.extend(configs)
                completed += 1
                elapsed = time.perf_counter() - started
                projected = elapsed / completed * (total_jobs - completed)
                print(f"checkpoint={cp.name} runtime={checkpoint['job_runtime_seconds']} projected_remaining={projected}", flush=True)
                if projected > a.max_projected_remaining_seconds:
                    reason = f"Stopped because projected remaining runtime {projected} seconds exceeded the pre-fixed 10800-second threshold."
                    write_outputs(out, records, a, started, "stopped", reason, projected)
                    print(reason, flush=True)
                    return
                write_outputs(out, records, a, started, "running", projected=projected)
    write_outputs(out, records, a, started, "completed")


if __name__ == "__main__":
    main()
