#!/usr/bin/env python3
"""Fold-local IPCW metrics and cross-fitted calibration for Full/Clinical-only OOF risks."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc, integrated_brier_score
from sksurv.util import Surv

MODELS=("coxph","rsf","deepsurv"); H=np.array([365.,1095.,1825.]); GRID=np.linspace(365.,1825.,100)

def evaluate(d: pd.DataFrame) -> dict:
    output={}; survival_estimates={m:np.full((len(d),len(GRID)),np.nan) for m in MODELS}
    for m in MODELS:
        folds=[]
        for held in sorted(d.fold.unique()):
            train=d[d.fold.ne(held)]; test=d[d.fold.eq(held)]; st=Surv.from_arrays(train.event.astype(bool),train.time); se=Surv.from_arrays(test.event.astype(bool),test.time)
            uno=float(concordance_index_ipcw(st,se,test[m].to_numpy())[0]); auc,weighted=cumulative_dynamic_auc(st,se,test[m].to_numpy(),H)
            mu=float(train[m].mean());sd=float(train[m].std(ddof=0));sd=sd if sd>0 and np.isfinite(sd) else 1.0
            cal=CoxPHFitter(penalizer=.1); z=train[[m,"time","event"]].rename(columns={m:"risk"});z['risk']=(z['risk']-mu)/sd;cal.fit(z,"time","event")
            ztest=pd.DataFrame({'risk':(test[m]-mu)/sd});sf=cal.predict_survival_function(ztest,times=GRID).T.to_numpy();survival_estimates[m][test.index,:]=sf
            ibs=float(integrated_brier_score(st,se,sf,GRID))
            folds.append({"fold":int(held),"n":len(test),"uno_c_index":uno,"auc_1y":float(auc[0]),"auc_3y":float(auc[1]),"auc_5y":float(auc[2]),"weighted_mean_auc":float(weighted),"ibs":ibs})
        output[m]={"folds":folds,**{f"mean_{k}":float(np.average([x[k] for x in folds],weights=[x['n'] for x in folds])) for k in ("uno_c_index","auc_1y","auc_3y","auc_5y","weighted_mean_auc","ibs")}}
        pred=survival_estimates[m][:,-1]; bins=np.asarray(pd.qcut(pred,5,labels=False,duplicates="drop"));groups=[]
        for level in np.unique(bins[np.isfinite(bins)]):
            pos=np.flatnonzero(bins==level);g=d.iloc[pos];km=KaplanMeierFitter().fit(g.time.astype(float),event_observed=g.event.astype(int))
            s=float(km.predict(1825.));ci=km.confidence_interval_survival_function_;row=ci.loc[ci.index[ci.index<=1825.][-1]] if np.any(ci.index<=1825.) else pd.Series([1.,1.])
            groups.append({"predicted_survival_mean":float(pred[pos].mean()),"observed_km_survival":s,"greenwood_95_ci":[float(row.iloc[0]),float(row.iloc[1])],"n_records":len(g),"n_unique_participants":int(g.sample_key.str[:12].nunique()),"at_risk_unique_participants_at_5y":int(g.loc[g.time>=1825.,'sample_key'].str[:12].nunique()),"events_by_5y":int(((g.event==1)&(g.time<=1825.)).sum())})
        output[m]["calibration_5year_groups"]=groups
    return output

def main():
    here=Path(__file__).resolve().parent;repo=here.parents[1];target=repo/'results/primary_tumor_only';corrected=pd.read_parquet(repo/'results/phase7/phase7_deepsurv_oof.parquet');result={"method":{ "time_dependent_auc_horizons_days":H.tolist(),"mean_auc_definition":"sksurv cumulative_dynamic_auc survival-function-weighted mean over the supplied 1/3/5-year grid; not arithmetic mean", "ibs_grid":"100 equally spaced points from day 365 through day 1825", "ipcw_censoring_source":"for each held-out outer fold, the other four OOF-fold outcomes; never that held fold", "aggregation":"test-fold-N-weighted mean", "absolute_survival":"one-covariate Cox calibration fitted on the other four OOF folds", "coxph_rsf_source":"retained Phase 2/5 ablation OOF predictions", "deepsurv_source":"Phase 7 leakage-corrected OOF predictions from results/phase7/phase7_deepsurv_oof.parquet"}}
    paths={"mixed":repo/'results/mixed_cohort/ablation_oof_extended.parquet',"primary":repo/'results/primary_tumor_only/ablation_oof_primary.parquet'}
    for cohort,path in paths.items():
        all_=pd.read_parquet(path);all_=all_[all_.configuration.isin(('Clinical-only','Full'))].copy();repl=corrected[corrected.cohort.eq(cohort)][['sample_key','fold','configuration','deepsurv']].rename(columns={'deepsurv':'deepsurv_corrected'});all_=all_.merge(repl,on=['sample_key','fold','configuration'],how='left',validate='one_to_one')
        if all_.deepsurv_corrected.isna().any():raise RuntimeError(f'Missing corrected DeepSurv predictions for {cohort}')
        all_['deepsurv']=all_.pop('deepsurv_corrected');result[cohort]={}
        for config in ("Clinical-only","Full"):
            d=all_[all_.configuration.eq(config)].reset_index(drop=True);print(cohort,config,len(d));result[cohort][config]=evaluate(d)
    target.mkdir(parents=True,exist_ok=True);(target/'additional_metrics_v2.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    lines=['# Additional Metrics v2','',result['method']['mean_auc_definition']+'.',f"IBS grid: {result['method']['ibs_grid']}.",f"IPCW censoring source: {result['method']['ipcw_censoring_source']}.",'']
    for cohort,sets in result.items():
        if cohort=='method':continue
        lines += [f"## {cohort}",'','| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |','|---|---|---:|---:|---:|---:|---:|---:|']
        for config,models in sets.items():
            for m,x in models.items():lines.append(f"| {config} | {m} | {x['mean_uno_c_index']} | {x['mean_auc_1y']} | {x['mean_auc_3y']} | {x['mean_auc_5y']} | {x['mean_weighted_mean_auc']} | {x['mean_ibs']} |")
        lines += ['', 'Five-year calibration group details (including Greenwood CI, unique participants at risk, and events) are preserved without rounding in `additional_metrics_v2.json`.', '']
    (target/'additional_metrics_v2.md').write_text('\n'.join(lines),encoding='utf-8');print(f'additional_metrics_v2.md/json written to {target}')
if __name__=='__main__':main()
