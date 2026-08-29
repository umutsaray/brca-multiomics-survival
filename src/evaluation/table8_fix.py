#!/usr/bin/env python3
"""Recompute Table 8 with corrected DeepSurv and fold-stratified bootstrap."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored

MODELS = ("deepsurv", "coxph", "rsf")
COMPARISONS = {"deepsurv_vs_coxph": "coxph", "deepsurv_vs_rsf": "rsf"}


def cindex(event, time, risk):
    return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])


def load_current(repo: Path) -> pd.DataFrame:
    classical = pd.read_parquet(repo / "results/mixed_cohort/ablation_oof_extended.parquet")
    classical = classical[classical.configuration.eq("Full")][["sample_key", "fold", "event", "time", "coxph", "rsf"]]
    dl = pd.read_parquet(repo / "results/phase7/phase7_deepsurv_oof.parquet")
    dl = dl[(dl.cohort.eq("mixed")) & (dl.configuration.eq("Full"))][["sample_key", "fold", "event", "time", "deepsurv"]]
    data = classical.merge(dl, on=["sample_key", "fold", "event", "time"], validate="one_to_one")
    if len(data) != len(classical) or data[list(MODELS)].isna().any().any():
        raise RuntimeError("Incomplete alignment of corrected DeepSurv and classical Full OOF predictions")
    data["participant"] = data.sample_key.str[:12]
    return data


def prepare_folds(data: pd.DataFrame):
    prepared = {}
    for fold in sorted(data.fold.unique()):
        fd = data[data.fold.eq(fold)].sort_values("sample_key").reset_index(drop=True)
        groups = fd.groupby("participant", sort=True).indices
        participants = np.array(list(groups))
        prepared[int(fold)] = {
            "n_records": len(fd), "n_participants": len(participants),
            "groups": [np.asarray(groups[p], dtype=int) for p in participants],
            "event": fd.event.to_numpy(bool), "time": fd.time.to_numpy(float),
            **{m: fd[m].to_numpy(float) for m in MODELS},
        }
    return prepared


def observed(prepared):
    fold_rows=[]
    for fold,item in prepared.items():
        cs={m:cindex(item["event"],item["time"],item[m]) for m in MODELS}
        fold_rows.append({"fold":fold,"n_records":item["n_records"],"n_participants":item["n_participants"],
                          "cindex":cs,"deepsurv_minus_coxph":cs["deepsurv"]-cs["coxph"],
                          "deepsurv_minus_rsf":cs["deepsurv"]-cs["rsf"]})
    return fold_rows


def bootstrap(prepared, iterations, seed):
    rng=np.random.default_rng(seed); arrays={k:np.empty(iterations) for k in COMPARISONS}
    for iteration in range(iterations):
        deltas={k:[] for k in COMPARISONS}
        for item in prepared.values():
            draw=rng.integers(0,item["n_participants"],size=item["n_participants"])
            idx=np.concatenate([item["groups"][i] for i in draw])
            scores={m:cindex(item["event"][idx],item["time"][idx],item[m][idx]) for m in MODELS}
            for name,other in COMPARISONS.items(): deltas[name].append(scores["deepsurv"]-scores[other])
        for name in COMPARISONS: arrays[name][iteration]=np.mean(deltas[name])
        if (iteration+1)%1000==0: print(f"bootstrap {iteration+1}/{iterations}",flush=True)
    return arrays


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--iterations",type=int,default=9999);p.add_argument("--seed",type=int,default=42);a=p.parse_args()
    repo=Path(__file__).resolve().parents[2];out=repo/"results/phase8";out.mkdir(parents=True,exist_ok=True)
    data=load_current(repo);prepared=prepare_folds(data);folds=observed(prepared);arrays=bootstrap(prepared,a.iterations,a.seed)
    old=json.loads((repo/"results/mixed_cohort/statistical_comparison_v2.json").read_text(encoding="utf-8"))
    result={"cohort":"mixed","configuration":"Full","n_records":len(data),"n_participants":int(data.participant.nunique()),
            "method":"Harrell C within each outer fold; unweighted mean of five fold deltas; participant bootstrap independently within each fold",
            "iterations":a.iterations,"seed":a.seed,"folds":folds,"comparisons":{}}
    for name,other in COMPARISONS.items():
        obs=np.array([f[f"deepsurv_minus_{other}"] for f in folds]);vals=arrays[name];lo,hi=np.percentile(vals,[2.5,97.5]);p0=min(1.,2*min((np.count_nonzero(vals<=0)+1)/(a.iterations+1),(np.count_nonzero(vals>=0)+1)/(a.iterations+1)))
        result["comparisons"][name]={"mean_fold_delta":float(obs.mean()),"sample_sd_fold_delta":float(obs.std(ddof=1)),"ci_95":[float(lo),float(hi)],"two_sided_p":float(p0),"contains_zero":bool(lo<=0<=hi),"bootstrap_deltas":vals.tolist()}
    (out/"table8_fix_results.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    lines=["# Table 8 Phase 8 correction","",
           "Mixed cohort, Full multi-omics configuration. DeepSurv risks are the Phase 7 leakage-corrected out-of-fold predictions; CoxPH and RSF risks are the unchanged retained Full-configuration predictions.","",
           f"New method: Harrell C is calculated separately inside each outer fold. Fold delta is DeepSurv C minus comparator C. The point estimate is the unweighted mean of five fold deltas. In each of {a.iterations:,} iterations, participants are resampled with replacement independently inside every fold; all three models use identical resample indices. Two-sided p-values use the paired bootstrap distribution with a +1 finite-sample correction. Seed: {a.seed}.","",
           "## Old pooled versus current fold-stratified results","",
           "| Comparison | Old pooled/leaky delta | Old 95% CI | Old p | New fold-stratified/corrected delta | Fold SD | New 95% CI | New p | Direction changed? | Significance changed? |","|---|---:|---:|---:|---:|---:|---:|---:|---|---|"]
    for name,other in COMPARISONS.items():
        o=old["comparisons"][name];n=result["comparisons"][name];old_sig=not o["contains_zero"];new_sig=not n["contains_zero"];direction_changed=np.sign(o["point_difference"])!=np.sign(n["mean_fold_delta"])
        lines.append(f"| DeepSurv − {other.upper()} | {o['point_difference']:.6f} | [{o['ci_95'][0]:.6f}, {o['ci_95'][1]:.6f}] | {o['paired_bootstrap_p_value']:.6f} | {n['mean_fold_delta']:.6f} | {n['sample_sd_fold_delta']:.6f} | [{n['ci_95'][0]:.6f}, {n['ci_95'][1]:.6f}] | {n['two_sided_p']:.6f} | {'Yes' if direction_changed else 'No'} | {'Yes' if old_sig!=new_sig else 'No'} |")
    lines += ["","## Fold-level current values","","| Fold | DeepSurv C | CoxPH C | RSF C | DeepSurv − CoxPH | DeepSurv − RSF |","|---:|---:|---:|---:|---:|---:|"]
    for f in folds: lines.append(f"| {f['fold']} | {f['cindex']['deepsurv']:.6f} | {f['cindex']['coxph']:.6f} | {f['cindex']['rsf']:.6f} | {f['deepsurv_minus_coxph']:.6f} | {f['deepsurv_minus_rsf']:.6f} |")
    lines += ["", "## Observed change", ""]
    for name,other in COMPARISONS.items():
        o=old["comparisons"][name];n=result["comparisons"][name];old_sig=not o["contains_zero"];new_sig=not n["contains_zero"]
        direction="changed" if np.sign(o["point_difference"])!=np.sign(n["mean_fold_delta"]) else "did not change"
        significance="changed" if old_sig!=new_sig else "did not change"
        lines.append(f"- DeepSurv − {other.upper()}: direction {direction}; significance status {significance}. Old and new significance statuses are {'significant' if old_sig else 'not significant'} and {'significant' if new_sig else 'not significant'}, respectively.")
    (out/"table8_fix_results.md").write_text("\n".join(lines),encoding="utf-8");print(f"Wrote {out/'table8_fix_results.md'} and JSON companion")

if __name__=="__main__":main()
