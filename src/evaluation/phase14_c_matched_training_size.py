#!/usr/bin/env python3
"""Phase 14 Check C: refit CoxPH/RSF on the exact Phase-7 85% fit subsets."""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import lifelines
import numpy as np
import pandas as pd
import sklearn
import sksurv
from lifelines import CoxPHFitter
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

HERE = Path(__file__).resolve().parent
for sub in ("common", "data_prep", "preprocessing", "feature_selection", "models", "evaluation", "explainability"):
    sys.path.insert(0, str(HERE.parents[0] / sub))
from phase1_common import grouped_folds, load_cohort  # noqa: E402
from phase7_deepsurv import inner_split  # noqa: E402


CONFIGS = ("Clinical-only", "Full")
MODELS = ("coxph", "rsf")


def cindex(event, duration, risk) -> float:
    return float(concordance_index_censored(np.asarray(event, bool), np.asarray(duration, float), np.asarray(risk, float))[0])


def fit_classical(xfit, yfit, xtest, ytest, seed, cox_penalizer, cox_max_steps, rsf_trees):
    out = {}
    event = yfit.event.astype(bool)
    event_var = xfit.loc[event].var()
    censor_var = xfit.loc[~event].var()
    cox_cols = [c for c in xfit if event_var[c] > 1e-12 and censor_var[c] > 1e-12]
    scaler = StandardScaler()
    a = scaler.fit_transform(xfit[cox_cols])
    b = scaler.transform(xtest[cox_cols])
    design = pd.DataFrame(a, columns=cox_cols)
    design["time"] = yfit.time.to_numpy()
    design["event"] = yfit.event.to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = CoxPHFitter(penalizer=cox_penalizer)
        model.fit(design, duration_col="time", event_col="event", fit_options={"max_steps": cox_max_steps})
    cox_risk = model.predict_partial_hazard(pd.DataFrame(b, columns=cox_cols)).to_numpy().ravel()
    out["coxph"] = {
        "risk": cox_risk,
        "cindex": cindex(ytest.event, ytest.time, cox_risk),
        "n_features_input": int(xfit.shape[1]),
        "n_features_after_event_censor_variance_filter": int(len(cox_cols)),
        "warnings": [str(w.message) for w in caught],
    }
    rsf = RandomSurvivalForest(
        n_estimators=rsf_trees,
        min_samples_split=10,
        min_samples_leaf=5,
        random_state=seed,
        n_jobs=-1,
    )
    rsf.fit(xfit.to_numpy("float32"), Surv.from_arrays(yfit.event.astype(bool), yfit.time))
    rsf_risk = rsf.predict(xtest.to_numpy("float32"))
    out["rsf"] = {
        "risk": rsf_risk,
        "cindex": cindex(ytest.event, ytest.time, rsf_risk),
        "n_features_input": int(xfit.shape[1]),
        "n_features_after_event_censor_variance_filter": None,
        "warnings": [],
    }
    return out


def summarize(values):
    array = np.asarray(values, float)
    return {"fold_cindices": array.tolist(), "mean": float(array.mean()), "sample_sd": float(array.std(ddof=1))}


def bootstrap_difference(data, config, model, iterations, seed):
    rng = np.random.default_rng(seed)
    folds = sorted(data.fold.unique())
    prepared = {}
    for fold in folds:
        d = data[(data.fold == fold) & (data.configuration == config)].sort_values("sample_key").reset_index(drop=True)
        d["participant"] = d.sample_key.str[:12]
        groups = d.groupby("participant", sort=True).indices
        participants = np.asarray(list(groups))
        prepared[fold] = {
            "groups": [np.asarray(groups[p], int) for p in participants],
            "event": d.event.to_numpy(bool),
            "time": d.time.to_numpy(float),
            "matched": d[f"matched_{model}"].to_numpy(float),
            "full": d[f"full_{model}"].to_numpy(float),
        }
    draws = np.empty(iterations, float)
    for i in range(iterations):
        fold_delta = []
        for fold in folds:
            item = prepared[fold]
            n = len(item["groups"])
            draw = rng.integers(0, n, size=n)
            idx = np.concatenate([item["groups"][j] for j in draw])
            matched = cindex(item["event"][idx], item["time"][idx], item["matched"][idx])
            full = cindex(item["event"][idx], item["time"][idx], item["full"][idx])
            fold_delta.append(matched - full)
        draws[i] = np.mean(fold_delta)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    p_lo = (np.count_nonzero(draws <= 0) + 1) / (iterations + 1)
    p_hi = (np.count_nonzero(draws >= 0) + 1) / (iterations + 1)
    return {
        "estimand": "unweighted mean of five fold-specific (matched-85%-fit C-index minus full-outer-training C-index)",
        "observed_mean_difference": float(np.mean([
            cindex(g.event, g.time, g[f"matched_{model}"]) - cindex(g.event, g.time, g[f"full_{model}"])
            for _, g in data[data.configuration == config].groupby("fold", sort=True)
        ])),
        "ci_95": [float(lo), float(hi)],
        "two_sided_p": float(min(1.0, 2 * min(p_lo, p_hi))),
        "iterations": iterations,
        "seed": seed,
        "bootstrap_differences": draws.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--repo-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--outer-seed", type=int, default=42)
    parser.add_argument("--bootstrap-iterations", type=int, default=9999)
    parser.add_argument("--bootstrap-seed", type=int, default=42)
    parser.add_argument("--cox-penalizer", type=float, default=0.1)
    parser.add_argument("--cox-max-steps", type=int, default=50)
    parser.add_argument("--rsf-trees", type=int, default=500)
    args = parser.parse_args()
    started = time.perf_counter()
    root, repo, out = args.input_dir.resolve(), args.repo_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    _, y = load_cohort(root, "primary")
    existing = json.loads((repo / "results/phase9/ablation_results_phase9.json").read_text(encoding="utf-8"))["primary"]
    current_oof = pd.read_parquet(repo / "results/phase9/ablation_oof_phase9.parquet")
    current_oof = current_oof[(current_oof.cohort == "primary") & current_oof.configuration.isin(CONFIGS)].copy()
    current_oof = current_oof[["sample_key", "fold", "configuration", "event", "time", "coxph", "rsf", "deepsurv"]]
    current_oof = current_oof.rename(columns={"coxph": "full_coxph", "rsf": "full_rsf"})

    matched_rows, fold_records = [], []
    matched_values = {config: {model: [] for model in MODELS} for config in CONFIGS}
    for fold, (outer_train, outer_test) in enumerate(grouped_folds(y, 5, args.outer_seed), 1):
        internal_seed = 4200 + fold
        fit_idx, validation_idx = inner_split(y, outer_train, internal_seed)
        fit_keys, test_keys = y.index[fit_idx], y.index[outer_test]
        stem = f"primary_fold{fold}_innerseed{internal_seed}"
        cache = root / "phase7_inner_cache"
        fit_path, test_path, meta_path = cache / f"{stem}_fit.parquet", cache / f"{stem}_test.parquet", cache / f"{stem}_meta.json"
        missing = [str(p) for p in (fit_path, test_path, meta_path) if not p.is_file()]
        if missing:
            raise FileNotFoundError("Required leakage-corrected Phase 7 cache missing: " + ", ".join(missing))
        xfit, xtest = pd.read_parquet(fit_path), pd.read_parquet(test_path)
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if not xfit.index.equals(fit_keys) or not xtest.index.equals(test_keys):
            raise RuntimeError(f"Phase 7 cache/index mismatch in fold {fold}")
        yfit, ytest = y.loc[fit_keys], y.loc[test_keys]
        record = {
            "fold": fold,
            "outer_train_n": int(len(outer_train)),
            "fit_subset_n": int(len(fit_idx)),
            "internal_validation_n": int(len(validation_idx)),
            "outer_test_n": int(len(outer_test)),
            "outer_seed": args.outer_seed,
            "internal_split_seed": internal_seed,
            "model_seed": args.outer_seed + fold,
            "phase7_cache_stem": stem,
            "selection": meta["selection"],
            "preprocessing": meta["preprocessing"],
            "configurations": {},
        }
        for config in CONFIGS:
            cols = [c for c in xfit if c.startswith("clin_")] if config == "Clinical-only" else list(xfit.columns)
            fits = fit_classical(
                xfit[cols], yfit, xtest[cols], ytest, args.outer_seed + fold,
                args.cox_penalizer, args.cox_max_steps, args.rsf_trees,
            )
            record["configurations"][config] = {}
            for model in MODELS:
                rec = {k: v for k, v in fits[model].items() if k != "risk"}
                record["configurations"][config][model] = rec
                matched_values[config][model].append(rec["cindex"])
            for j, key in enumerate(test_keys):
                matched_rows.append({
                    "sample_key": key, "fold": fold, "configuration": config,
                    "event": int(ytest.iloc[j].event), "time": float(ytest.iloc[j].time),
                    "matched_coxph": float(fits["coxph"]["risk"][j]),
                    "matched_rsf": float(fits["rsf"]["risk"][j]),
                })
            print(f"fold={fold} config={config} coxph={fits['coxph']['cindex']} rsf={fits['rsf']['cindex']}", flush=True)
        fold_records.append(record)

    matched = pd.DataFrame(matched_rows)
    comparison = current_oof.merge(
        matched,
        on=["sample_key", "fold", "configuration", "event", "time"],
        how="inner",
        validate="one_to_one",
    )
    if len(comparison) != 2 * len(y):
        raise RuntimeError(f"OOF merge produced {len(comparison)} rows, expected {2 * len(y)}")

    results = {}
    for config in CONFIGS:
        results[config] = {}
        for model in MODELS:
            matched_summary = summarize(matched_values[config][model])
            full_summary = existing[config][model]
            boot = bootstrap_difference(comparison, config, model, args.bootstrap_iterations, args.bootstrap_seed)
            results[config][model] = {"matched_85_percent": matched_summary, "full_outer_training_table7": full_summary, "matched_minus_full_bootstrap": boot}

    deep = {config: existing[config]["deepsurv"] for config in CONFIGS}
    ordering = {"clinical_vs_full": {}, "deepsurv_vs_classical": {}}
    for model in MODELS:
        before_delta = existing["Clinical-only"][model]["mean"] - existing["Full"][model]["mean"]
        after_delta = results["Clinical-only"][model]["matched_85_percent"]["mean"] - results["Full"][model]["matched_85_percent"]["mean"]
        ordering["clinical_vs_full"][model] = {
            "full_training_clinical_minus_full": before_delta,
            "matched_training_clinical_minus_full": after_delta,
            "ordering_changed": bool(np.sign(before_delta) != np.sign(after_delta)),
        }
    for config in CONFIGS:
        ordering["deepsurv_vs_classical"][config] = {}
        for model in MODELS:
            before_delta = deep[config]["mean"] - existing[config][model]["mean"]
            after_delta = deep[config]["mean"] - results[config][model]["matched_85_percent"]["mean"]
            ordering["deepsurv_vs_classical"][config][model] = {
                "deepsurv_mean": deep[config]["mean"],
                "classical_full_training_mean": existing[config][model]["mean"],
                "classical_matched_training_mean": results[config][model]["matched_85_percent"]["mean"],
                "before_deepsurv_minus_classical": before_delta,
                "after_deepsurv_minus_classical": after_delta,
                "ordering_changed": bool(np.sign(before_delta) != np.sign(after_delta)),
            }

    payload = {
        "check": "C - matched training-set size for CoxPH / RSF vs DeepSurv",
        "status": "completed",
        "method": {
            "cohort": "corrected primary-tumor-only",
            "n_participants": int(len(y)),
            "n_events": int(y.event.sum()),
            "outer_split": "5-fold StratifiedGroupKFold, participant-grouped/event-stratified",
            "outer_seed": args.outer_seed,
            "internal_split": "same Phase 7 first held fold of 7-fold StratifiedGroupKFold; fit subset only",
            "internal_seed_schema": "4200 + outer_fold",
            "feature_selection_imputation_scaling": "read from Phase 7 cache fitted only on internal-train; Cox scaling refit on same internal-train",
            "cox_penalizer": args.cox_penalizer,
            "cox_max_steps": args.cox_max_steps,
            "rsf_trees": args.rsf_trees,
            "rsf_min_samples_split": 10,
            "rsf_min_samples_leaf": 5,
            "model_seed_schema": "42 + outer_fold",
            "bootstrap": "participant resampling within each outer fold; unweighted mean of five fold-specific differences",
            "bootstrap_iterations": args.bootstrap_iterations,
            "bootstrap_seed": args.bootstrap_seed,
            "deviation_from_requested_wording": "Phase 7 uses one of seven grouped folds for internal validation (~14.3%), not an independently generated exact 15% split; the exact existing Phase 7 split and cache were reused.",
        },
        "command": " ".join([sys.executable, *sys.argv]),
        "versions": {
            "python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
            "scikit_learn": sklearn.__version__, "scikit_survival": sksurv.__version__, "lifelines": lifelines.__version__,
        },
        "runtime_seconds": None,
        "folds": fold_records,
        "results": results,
        "ordering_checks": ordering,
    }
    payload["runtime_seconds"] = time.perf_counter() - started
    (out / "phaseC_matched_training_size.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# Phase 14 Check C — Matched training-set size",
        "",
        "Status: **completed**.",
        "",
        "The exact Phase 7 internal fit subsets and cached preprocessing/feature selection were reused. Phase 7 holds out one of seven grouped folds (~14.3%), so the retained fit subset is described as matched-85% but is not a newly generated exact 15% holdout.",
        "",
        "## C-index results",
        "",
        "| Configuration | Model | Full outer-training mean | Full SD | Matched fit-subset mean | Matched SD | Matched fold C-index values |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for config in CONFIGS:
        for model in MODELS:
            r = results[config][model]
            md.append(
                f"| {config} | {model} | {r['full_outer_training_table7']['mean']} | {r['full_outer_training_table7']['sample_sd']} | "
                f"{r['matched_85_percent']['mean']} | {r['matched_85_percent']['sample_sd']} | `{json.dumps(r['matched_85_percent']['fold_cindices'])}` |"
            )
    md += ["", "## Fold-stratified paired bootstrap: matched minus full", "", "| Configuration | Model | Mean difference | 95% CI | Two-sided p | Iterations | Seed |", "|---|---|---:|---|---:|---:|---:|"]
    for config in CONFIGS:
        for model in MODELS:
            b = results[config][model]["matched_minus_full_bootstrap"]
            md.append(f"| {config} | {model} | {b['observed_mean_difference']} | `{json.dumps(b['ci_95'])}` | {b['two_sided_p']} | {b['iterations']} | {b['seed']} |")
    md += ["", "## Ordering checks", "", "Clinical-only versus Full:", ""]
    for model, rec in ordering["clinical_vs_full"].items():
        md.append(f"- `{model}`: changed = `{str(rec['ordering_changed']).lower()}`; full-training Clinical−Full = `{rec['full_training_clinical_minus_full']}`; matched Clinical−Full = `{rec['matched_training_clinical_minus_full']}`.")
    md += ["", "DeepSurv versus classical model:", ""]
    for config, models in ordering["deepsurv_vs_classical"].items():
        for model, rec in models.items():
            md.append(f"- `{config}`, DeepSurv vs `{model}`: changed = `{str(rec['ordering_changed']).lower()}`; before DeepSurv−classical = `{rec['before_deepsurv_minus_classical']}`; after = `{rec['after_deepsurv_minus_classical']}`.")
    md += ["", "## Fold details", ""]
    for fold in fold_records:
        md.append(f"- Fold `{fold['fold']}`: outer train `{fold['outer_train_n']}`, fit subset `{fold['fit_subset_n']}`, internal validation `{fold['internal_validation_n']}`, outer test `{fold['outer_test_n']}`, internal seed `{fold['internal_split_seed']}`, model seed `{fold['model_seed']}`.")
    md += [
        "", "## Reproducibility", "",
        f"Command: `{payload['command']}`  ",
        f"Versions: `{json.dumps(payload['versions'])}`  ",
        f"Runtime seconds: `{payload['runtime_seconds']}`.",
    ]
    (out / "phaseC_matched_training_size.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "runtime_seconds": payload["runtime_seconds"], "results": {c: {m: results[c][m]["matched_85_percent"] for m in MODELS} for c in CONFIGS}, "ordering": ordering}, indent=2))


if __name__ == "__main__":
    main()
