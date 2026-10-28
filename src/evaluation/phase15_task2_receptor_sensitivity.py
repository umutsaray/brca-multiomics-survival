#!/usr/bin/env python3
"""Phase 15 Task 2: receptor-augmented survival sensitivity analyses."""
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
import statsmodels
import torch
from lifelines import CoxPHFitter
from sklearn.compose import ColumnTransformer
from sklearn.impute import KNNImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
for sub in ("common", "data_prep", "preprocessing", "feature_selection", "models", "evaluation", "explainability"):
    sys.path.insert(0, str(HERE.parents[0] / sub))
from impute_and_split import scan_train_methylation_variance  # noqa: E402
from phase1_common import (  # noqa: E402
    fold_clinical, grouped_folds, load_cohort, load_selected_methylation,
    prepare_fold, rank_fold, select_from_rankings,
)
from phase7_deepsurv import fit_network, inner_split, prepare_inner  # noqa: E402


CONFIGS = ("Clinical-only", "Clinical + ER/PR/HER2", "Full", "Full + ER/PR/HER2")
MODELS = ("coxph", "rsf", "deepsurv")
COMPARISONS = (
    ("Clinical + ER/PR/HER2 minus Clinical-only", "Clinical + ER/PR/HER2", "Clinical-only"),
    ("Clinical-only minus Full", "Clinical-only", "Full"),
    ("Clinical + ER/PR/HER2 minus Full", "Clinical + ER/PR/HER2", "Full"),
    ("Full + ER/PR/HER2 minus Full", "Full + ER/PR/HER2", "Full"),
)
RECEPTORS = ("ER", "PR", "HER2")


def cindex(event, duration, risk) -> float:
    return float(concordance_index_censored(np.asarray(event, bool), np.asarray(duration, float), np.asarray(risk, float))[0])


def receptor_frame(task1: dict, sample_keys: pd.Index) -> pd.DataFrame:
    mapping = task1["participant_model_categories"]
    rows = [{r: mapping[key[:12]][r] for r in RECEPTORS} for key in sample_keys]
    return pd.DataFrame(rows, index=sample_keys)


def encode_receptors(task1: dict, fit_keys: pd.Index, *other_keys: pd.Index):
    all_keys = fit_keys.append(list(other_keys)) if other_keys else fit_keys
    raw = receptor_frame(task1, all_keys)
    encoder = OneHotEncoder(
        categories=[["Negative", "Positive", "Not available"]] * 3,
        drop="first", handle_unknown="ignore", sparse_output=False, dtype=np.float32,
    )
    encoder.fit(raw.loc[fit_keys, list(RECEPTORS)])
    names = ["clin_receptor_" + name for name in encoder.get_feature_names_out(RECEPTORS)]
    def transform(keys):
        return pd.DataFrame(encoder.transform(raw.loc[keys, list(RECEPTORS)]), index=keys, columns=names, dtype="float32")
    return (transform(fit_keys), *(transform(keys) for keys in other_keys)), names


def matrix_diagnostics(matrix: np.ndarray) -> dict:
    rank = int(np.linalg.matrix_rank(matrix))
    return {
        "n_rows": int(matrix.shape[0]), "n_columns": int(matrix.shape[1]),
        "matrix_rank": rank, "rank_deficient": bool(rank < matrix.shape[1]),
    }


def fit_classical(xfit, yfit, xtest, ytest, seed, args):
    output = {}
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
        cox = CoxPHFitter(penalizer=args.cox_penalizer)
        cox.fit(design, duration_col="time", event_col="event", fit_options={"max_steps": args.cox_max_steps})
    cox_risk = cox.predict_partial_hazard(pd.DataFrame(b, columns=cox_cols)).to_numpy().ravel()
    output["coxph"] = {
        "risk": cox_risk, "cindex": cindex(ytest.event, ytest.time, cox_risk),
        "diagnostics": {**matrix_diagnostics(a), "input_columns": int(xfit.shape[1]), "columns_after_variance_filter": len(cox_cols), "convergence_warnings": [str(w.message) for w in caught]},
    }
    raw = xfit.to_numpy("float32")
    rsf = RandomSurvivalForest(
        n_estimators=args.rsf_trees, min_samples_split=10, min_samples_leaf=5,
        random_state=seed, n_jobs=-1,
    )
    rsf.fit(raw, Surv.from_arrays(yfit.event.astype(bool), yfit.time))
    rsf_risk = rsf.predict(xtest.to_numpy("float32"))
    output["rsf"] = {
        "risk": rsf_risk, "cindex": cindex(ytest.event, ytest.time, rsf_risk),
        "diagnostics": {**matrix_diagnostics(raw), "input_columns": int(xfit.shape[1]), "columns_after_variance_filter": None, "convergence_warnings": []},
    }
    return output


def prepare_outer_uncached(root, base, y, train_idx, test_idx, args):
    train_keys, test_keys = y.index[train_idx], y.index[test_idx]
    chosen, total, seconds = scan_train_methylation_variance(
        root / "TCGA-BRCA.methylation450.tsv.gz", train_keys, y.index,
        args.meth_prefilter, args.methylation_chunksize,
    )
    meth = load_selected_methylation(root / "TCGA-BRCA.methylation450.tsv.gz", chosen, y.index)
    imputer = KNNImputer(n_neighbors=args.knn_neighbors, keep_empty_features=True)
    mtrain = pd.DataFrame(imputer.fit_transform(meth.loc[train_keys]), index=train_keys, columns=chosen, dtype="float32")
    mtest = pd.DataFrame(imputer.transform(meth.loc[test_keys]), index=test_keys, columns=chosen, dtype="float32")
    ctrain, ctest = fold_clinical(root, y.index, train_keys, test_keys)
    return (
        pd.concat([base.loc[train_keys], mtrain, ctrain], axis=1),
        pd.concat([base.loc[test_keys], mtest, ctest], axis=1),
        {"cache": False, "methylation_total": total, "methylation_prefiltered": len(chosen), "variance_seconds": seconds, "knn_neighbors": args.knn_neighbors},
    )


def config_columns(frame: pd.DataFrame, receptor_names: list[str]) -> dict[str, list[str]]:
    receptors = set(receptor_names)
    clinical = [c for c in frame if c.startswith("clin_") and c not in receptors]
    base_full = [c for c in frame if c not in receptors]
    return {
        "Clinical-only": clinical,
        "Clinical + ER/PR/HER2": clinical + receptor_names,
        "Full": base_full,
        "Full + ER/PR/HER2": base_full + receptor_names,
    }


def summarize_fold_records(folds: list[dict]) -> dict:
    result = {}
    for config in CONFIGS:
        result[config] = {}
        for model in MODELS:
            values = [f["configurations"][config][model]["cindex"] for f in folds]
            result[config][model] = {
                "fold_cindices": values,
                "mean": float(np.mean(values)),
                "sample_sd": float(np.std(values, ddof=1)),
            }
    return result


def bootstrap_family(oof: pd.DataFrame, iterations: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    folds = sorted(oof.fold.unique())
    prepared = {}
    for fold in folds:
        by_config = {}
        for config in CONFIGS:
            d = oof[(oof.fold == fold) & (oof.configuration == config)].sort_values("sample_key").reset_index(drop=True)
            d["participant"] = d.sample_key.str[:12]
            groups = d.groupby("participant", sort=True).indices
            participants = np.asarray(list(groups))
            by_config[config] = {
                "participants": participants,
                "groups": [np.asarray(groups[p], int) for p in participants],
                "event": d.event.to_numpy(bool), "time": d.time.to_numpy(float),
                **{model: d[model].to_numpy(float) for model in MODELS},
            }
        reference = by_config[CONFIGS[0]]["participants"]
        if any(not np.array_equal(reference, by_config[c]["participants"]) for c in CONFIGS[1:]):
            raise RuntimeError(f"Participant mismatch across configurations in fold {fold}")
        prepared[fold] = by_config
    arrays = {(name, model): np.empty(iterations, float) for name, _, _ in COMPARISONS for model in MODELS}
    for i in range(iterations):
        fold_deltas = {key: [] for key in arrays}
        for fold in folds:
            by_config = prepared[fold]
            n = len(by_config[CONFIGS[0]]["participants"])
            draw = rng.integers(0, n, size=n)
            scores = {}
            for config in CONFIGS:
                item = by_config[config]
                idx = np.concatenate([item["groups"][j] for j in draw])
                for model in MODELS:
                    scores[(config, model)] = cindex(item["event"][idx], item["time"][idx], item[model][idx])
            for name, left, right in COMPARISONS:
                for model in MODELS:
                    fold_deltas[(name, model)].append(scores[(left, model)] - scores[(right, model)])
        for key in arrays:
            arrays[key][i] = np.mean(fold_deltas[key])
        if (i + 1) % 1000 == 0:
            print(f"bootstrap {i + 1}/{iterations}", flush=True)
    observed = {}
    for name, left, right in COMPARISONS:
        observed[name] = {}
        for model in MODELS:
            deltas = []
            for fold in folds:
                ld = oof[(oof.fold == fold) & (oof.configuration == left)]
                rd = oof[(oof.fold == fold) & (oof.configuration == right)]
                deltas.append(cindex(ld.event, ld.time, ld[model]) - cindex(rd.event, rd.time, rd[model]))
            values = arrays[(name, model)]
            lo, hi = np.percentile(values, [2.5, 97.5])
            p_lo = (np.count_nonzero(values <= 0) + 1) / (iterations + 1)
            p_hi = (np.count_nonzero(values >= 0) + 1) / (iterations + 1)
            observed[name][model] = {
                "fold_differences": deltas,
                "mean_difference": float(np.mean(deltas)),
                "ci_95": [float(lo), float(hi)],
                "raw_two_sided_p": float(min(1.0, 2 * min(p_lo, p_hi))),
                "bootstrap_differences": values.tolist(),
            }
    keys = [(name, model) for name, _, _ in COMPARISONS for model in MODELS]
    pvalues = [observed[name][model]["raw_two_sided_p"] for name, model in keys]
    adjusted = multipletests(pvalues, method="fdr_bh")[1]
    for (name, model), q in zip(keys, adjusted):
        observed[name][model]["bh_fdr_q_within_analysis_family"] = float(q)
    return observed


def make_markdown(title: str, payload: dict) -> str:
    lines = [
        f"# {title}", "",
        "All decision rules were fixed before model results were inspected. PAM50 was not used as a model input. This is an exploratory sensitivity analysis.", "",
        "## C-index", "",
        "| Configuration | Model | Mean | Sample SD | Fold values |", "|---|---|---:|---:|---|",
    ]
    for config in CONFIGS:
        for model in MODELS:
            rec = payload["summary"][config][model]
            lines.append(f"| {config} | {model} | {rec['mean']} | {rec['sample_sd']} | `{json.dumps(rec['fold_cindices'])}` |")
    lines += ["", "## Fold-stratified paired bootstrap", "", "| Comparison | Model | Mean difference | 95% CI | Raw p | BH-FDR q |", "|---|---|---:|---|---:|---:|"]
    for name, _, _ in COMPARISONS:
        for model in MODELS:
            rec = payload["bootstrap_comparisons"][name][model]
            lines.append(f"| {name} | {model} | {rec['mean_difference']} | `{json.dumps(rec['ci_95'])}` | {rec['raw_two_sided_p']} | {rec['bh_fdr_q_within_analysis_family']} |")
    lines += ["", "## Design-matrix and convergence diagnostics", "", "| Fold | Configuration | Model | Rows | Columns | Rank | Rank deficient | Convergence warnings |", "|---:|---|---|---:|---:|---:|---|---|"]
    for fold in payload["folds"]:
        for config in CONFIGS:
            for model in MODELS:
                d = fold["configurations"][config][model]["diagnostics"]
                lines.append(f"| {fold['fold']} | {config} | {model} | {d['n_rows']} | {d['n_columns']} | {d['matrix_rank']} | {str(d['rank_deficient']).lower()} | `{json.dumps(d['convergence_warnings'])}` |")
    if "event_counts_per_outer_test_fold" in payload:
        lines += ["", "## Complete-case fold event counts", "", f"Counts: `{json.dumps(payload['event_counts_per_outer_test_fold'])}`.", f"Folds with fewer than five test events: `{json.dumps(payload['folds_with_fewer_than_5_test_events'])}`."]
    if "reproduction_check_against_table7" in payload:
        lines += ["", "## Table 7 reproduction check", "", f"`{json.dumps(payload['reproduction_check_against_table7'])}`"]
    lines += ["", "## Reproducibility", "", f"Command: `{payload['command']}`  ", f"Versions: `{json.dumps(payload['versions'])}`  ", f"Seeds: `{json.dumps(payload['seeds'])}`  ", f"Runtime seconds: `{payload['runtime_seconds']}`."]
    return "\n".join(lines) + "\n"


def run_analysis(name, root, repo, task1, base, y, folds, args, use_main_cache, existing_oof=None):
    started = time.perf_counter()
    fold_records, rows = [], []
    for fold, (outer_train, outer_test) in enumerate(folds, 1):
        ytrain, ytest = y.iloc[outer_train], y.iloc[outer_test]
        if use_main_cache:
            xtrain, xtest, prep = prepare_fold(root, base, y, outer_train, outer_test, args.meth_prefilter, args.knn_neighbors, args.methylation_chunksize)
        else:
            xtrain, xtest, prep = prepare_outer_uncached(root, base, y, outer_train, outer_test, args)
        rankings = rank_fold(xtrain, ytrain)
        selected, selection = select_from_rankings(rankings, xtrain, args.top_n, True)
        selected = [c for c in selected if c in xtrain]
        xtrain, xtest = xtrain[selected], xtest[selected]
        (rtrain, rtest), receptor_names = encode_receptors(task1, ytrain.index, ytest.index)
        xtrain = pd.concat([xtrain, rtrain], axis=1)
        xtest = pd.concat([xtest, rtest], axis=1)
        configs = config_columns(xtrain, receptor_names)
        fr = {"fold": fold, "outer_train_n": len(outer_train), "outer_test_n": len(outer_test), "outer_test_events": int(ytest.event.sum()), "preprocessing": prep, "selection": selection, "configurations": {}}
        for config in CONFIGS:
            fits = fit_classical(xtrain[configs[config]], ytrain, xtest[configs[config]], ytest, args.outer_seed + fold, args)
            fr["configurations"][config] = {}
            for model in ("coxph", "rsf"):
                fr["configurations"][config][model] = {k: v for k, v in fits[model].items() if k != "risk"}
            for j, key in enumerate(ytest.index):
                rows.append({"sample_key": key, "fold": fold, "configuration": config, "event": int(ytest.iloc[j].event), "time": float(ytest.iloc[j].time), "coxph": float(fits['coxph']['risk'][j]), "rsf": float(fits['rsf']['risk'][j]), "deepsurv": np.nan})

        # Leakage-corrected DeepSurv: preprocessing, selection and receptor encoding fit on internal-train only.
        internal_seed = 4200 + fold
        fit_idx, val_idx = inner_split(y, outer_train, internal_seed)
        yfit, yval = y.iloc[fit_idx], y.iloc[val_idx]
        if use_main_cache:
            stem = f"primary_fold{fold}_innerseed{internal_seed}"
            cache = root / "phase7_inner_cache"
            xfit = pd.read_parquet(cache / f"{stem}_fit.parquet")
            xval = pd.read_parquet(cache / f"{stem}_validation.parquet")
            xdeep_test = pd.read_parquet(cache / f"{stem}_test.parquet")
            deep_meta = json.loads((cache / f"{stem}_meta.json").read_text(encoding="utf-8"))
        else:
            xfit, xval, xdeep_test, deep_prep = prepare_inner(root, base, y, fit_idx, val_idx, outer_test, args.meth_prefilter, args.knn_neighbors, args.methylation_chunksize)
            deep_rank = rank_fold(xfit, yfit)
            deep_selected, deep_selection = select_from_rankings(deep_rank, xfit, args.top_n, True)
            deep_selected = [c for c in deep_selected if c in xfit]
            xfit, xval, xdeep_test = xfit[deep_selected], xval[deep_selected], xdeep_test[deep_selected]
            deep_meta = {"selection": deep_selection, "preprocessing": deep_prep}
        (rf, rv, rt), deep_receptor_names = encode_receptors(task1, yfit.index, yval.index, ytest.index)
        xfit, xval, xdeep_test = pd.concat([xfit, rf], axis=1), pd.concat([xval, rv], axis=1), pd.concat([xdeep_test, rt], axis=1)
        deep_configs = config_columns(xfit, deep_receptor_names)
        for config in CONFIGS:
            if use_main_cache and config in {"Clinical-only", "Full"}:
                retained = existing_oof[(existing_oof.fold == fold) & (existing_oof.configuration == config)].set_index("sample_key").loc[ytest.index]
                risk = retained.deepsurv.to_numpy(float)
                score = cindex(ytest.event, ytest.time, risk)
                deep_info = {"cindex": score, "source": "retained current Phase 9 OOF prediction; not refitted", "best_epoch": None, "validation_cindex": None}
            else:
                fit = fit_network(xfit[deep_configs[config]], yfit, xval[deep_configs[config]], yval, xdeep_test[deep_configs[config]], ytest, args.outer_seed + fold, args.dl_epochs, args.dl_patience, args.dl_learning_rate, args.dl_weight_decay, args.dl_dropout)
                risk, score = fit["risk"], fit["cindex"]
                deep_info = {"cindex": score, "source": "Phase 9 leakage-corrected internal-validation code path", "best_epoch": fit["best_epoch"], "validation_cindex": fit["validation_cindex"]}
            matrix = xfit[deep_configs[config]].to_numpy("float32")
            deep_info["diagnostics"] = {**matrix_diagnostics(matrix), "input_columns": matrix.shape[1], "columns_after_variance_filter": None, "convergence_warnings": []}
            fr["configurations"][config]["deepsurv"] = deep_info
            index_map = {(row["sample_key"], row["fold"], row["configuration"]): row for row in rows[-4 * len(ytest):]}
            for j, key in enumerate(ytest.index):
                index_map[(key, fold, config)]["deepsurv"] = float(risk[j])
        fr["deepsurv_internal"] = {"internal_train_n": len(fit_idx), "internal_validation_n": len(val_idx), "internal_seed": internal_seed, "training_seed": args.outer_seed + fold, **deep_meta}
        fold_records.append(fr)
        print(name, "fold", fold, {c: {m: fr['configurations'][c][m]['cindex'] for m in MODELS} for c in CONFIGS}, flush=True)

    oof = pd.DataFrame(rows)
    if oof[list(MODELS)].isna().any().any():
        raise RuntimeError(f"{name}: missing OOF predictions")
    summary = summarize_fold_records(fold_records)
    boot = bootstrap_family(oof, args.bootstrap_iterations, args.bootstrap_seed)
    return fold_records, summary, boot, oof, time.perf_counter() - started


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--repo-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--outer-seed", type=int, default=42)
    parser.add_argument("--bootstrap-iterations", type=int, default=9999)
    parser.add_argument("--bootstrap-seed", type=int, default=42)
    parser.add_argument("--meth-prefilter", type=int, default=10000)
    parser.add_argument("--knn-neighbors", type=int, default=5)
    parser.add_argument("--methylation-chunksize", type=int, default=5000)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--cox-penalizer", type=float, default=0.1)
    parser.add_argument("--cox-max-steps", type=int, default=50)
    parser.add_argument("--rsf-trees", type=int, default=500)
    parser.add_argument("--dl-epochs", type=int, default=300)
    parser.add_argument("--dl-patience", type=int, default=30)
    parser.add_argument("--dl-learning-rate", type=float, default=0.001)
    parser.add_argument("--dl-weight-decay", type=float, default=0.0001)
    parser.add_argument("--dl-dropout", type=float, default=0.3)
    args = parser.parse_args()
    torch.set_num_threads(1)
    root, repo, out = args.input_dir.resolve(), args.repo_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    task1 = json.loads((out / "task1_receptor_missingness.json").read_text(encoding="utf-8"))
    base, y = load_cohort(root, "primary")
    retained_oof = pd.read_parquet(repo / "results/phase9/ablation_oof_phase9.parquet")
    retained_oof = retained_oof[(retained_oof.cohort == "primary") & retained_oof.configuration.isin(["Clinical-only", "Full"])]

    common = {
        "decision_rules": task1["decision_rules_fixed_before_model_results"],
        "exploratory": True,
        "pam50_used_as_model_input": False,
        "command": " ".join([sys.executable, *sys.argv]),
        "versions": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scikit_learn": sklearn.__version__, "scikit_survival": sksurv.__version__, "lifelines": lifelines.__version__, "statsmodels": statsmodels.__version__, "torch": torch.__version__},
        "seeds": {"outer": args.outer_seed, "deepsurv_internal_schema": "4200 + outer_fold", "model_schema": "42 + outer_fold", "bootstrap": args.bootstrap_seed},
        "bootstrap_iterations": args.bootstrap_iterations,
        "models": list(MODELS), "configurations": list(CONFIGS),
    }

    folds_a = grouped_folds(y, 5, args.outer_seed)
    fa, sa, ba, _, runtime_a = run_analysis("2A", root, repo, task1, base, y, folds_a, args, True, retained_oof)
    table7 = json.loads((repo / "results/phase9/ablation_results_phase9.json").read_text(encoding="utf-8"))["primary"]
    reproduction = {}
    for config in ("Clinical-only", "Full"):
        reproduction[config] = {}
        for model in ("coxph", "rsf", "deepsurv"):
            observed, expected = sa[config][model]["fold_cindices"], table7[config][model]["fold_cindices"]
            reproduction[config][model] = {"observed": observed, "table7": expected, "exact_match": observed == expected, "differences": [a - b for a, b in zip(observed, expected)]}
    payload_a = {**common, "task": "2A missing-indicator analysis", "status": "completed", "cohort_n": len(y), "events": int(y.event.sum()), "outer_folds": "identical to Table 7", "folds": fa, "summary": sa, "bootstrap_comparisons": ba, "reproduction_check_against_table7": reproduction, "runtime_seconds": runtime_a}
    (out / "task2a_missing_indicator.json").write_text(json.dumps(payload_a, indent=2), encoding="utf-8")
    (out / "task2a_missing_indicator.md").write_text(make_markdown("Phase 15 Task 2A — Missing-indicator receptor sensitivity", payload_a), encoding="utf-8")

    complete = set(task1["complete_case"]["participants"])
    mask = y.index.str[:12].isin(complete)
    base_b, y_b = base.loc[mask], y.loc[mask]
    folds_b = grouped_folds(y_b, 5, args.outer_seed)
    event_counts = [int(y_b.iloc[test].event.sum()) for _, test in folds_b]
    fb, sb, bb, _, runtime_b = run_analysis("2B", root, repo, task1, base_b, y_b, folds_b, args, False, None)
    payload_b = {**common, "task": "2B complete-case analysis", "status": "completed", "cohort_n": len(y_b), "events": int(y_b.event.sum()), "outer_folds": "new participant-grouped/event-stratified folds because the complete-case cohort differs from the main cohort", "event_counts_per_outer_test_fold": event_counts, "folds_with_fewer_than_5_test_events": [i + 1 for i, n in enumerate(event_counts) if n < 5], "folds": fb, "summary": sb, "bootstrap_comparisons": bb, "runtime_seconds": runtime_b}
    (out / "task2b_complete_case.json").write_text(json.dumps(payload_b, indent=2), encoding="utf-8")
    (out / "task2b_complete_case.md").write_text(make_markdown("Phase 15 Task 2B — Complete-case receptor sensitivity", payload_b), encoding="utf-8")
    print(json.dumps({"task2a_runtime_seconds": runtime_a, "task2b_runtime_seconds": runtime_b, "task2a_summary": sa, "task2b_summary": sb, "event_counts_2b": event_counts}, indent=2))


if __name__ == "__main__":
    main()
