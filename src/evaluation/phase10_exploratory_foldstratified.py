#!/usr/bin/env python3
"""Fold-stratified bootstrap for all 48 exploratory clinical-vs-omics tests."""
from __future__ import annotations
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from sksurv.metrics import concordance_index_censored
from statsmodels.stats.multitest import multipletests

MODELS=('coxph','rsf','deepsurv')
CONFIGS=('RNA-seq only','Methylation only','miRNA only','Clinical + RNA-seq','Clinical + RNA-seq + Methylation','Clinical + RNA-seq + miRNA','Clinical + Methylation','Clinical + miRNA')

def ci(e,t,r):return float(concordance_index_censored(e.astype(bool),t.astype(float),r.astype(float))[0])

def one_test(data,cohort,model,config,iterations,seed):
 prepared=[];fold_rows=[]
 for fold in sorted(data.fold.unique()):
  z=data[(data.fold.eq(fold)) & (data.configuration.eq(config))].dropna(subset=[model]).sort_values('sample_key').reset_index(drop=True)
  if z.empty:continue
  participants=z.sample_key.str[:12];groups=z.groupby(participants,sort=True).indices;blocks=[np.asarray(x,dtype=int) for x in groups.values()];e=z.event.to_numpy(bool);t=z.time.to_numpy(float);a=z[model].to_numpy(float);b=z[f'{model}_clinical'].to_numpy(float);ca=ci(e,t,a);cb=ci(e,t,b);fold_rows.append({'fold':int(fold),'n_records':len(z),'n_participants':len(blocks),'configuration_cindex':ca,'clinical_cindex':cb,'delta':ca-cb});prepared.append((blocks,e,t,a,b))
 if not prepared:raise RuntimeError(f'No estimable folds: {cohort} {model} {config}')
 rng=np.random.default_rng(seed);arr=np.empty(iterations)
 for i in range(iterations):
  ds=[]
  for blocks,e,t,a,b in prepared:
   draw=rng.choice(len(blocks),size=len(blocks),replace=True);idx=np.concatenate([blocks[j] for j in draw]);ds.append(ci(e[idx],t[idx],a[idx])-ci(e[idx],t[idx],b[idx]))
  arr[i]=np.mean(ds)
 obs=float(np.mean([x['delta'] for x in fold_rows]));lo,hi=np.percentile(arr,[2.5,97.5]);pv=min(1.,2*min((np.count_nonzero(arr<=0)+1)/(iterations+1),(np.count_nonzero(arr>=0)+1)/(iterations+1)))
 return {'cohort':cohort,'model':model,'configuration':config,'folds_used':[x['fold'] for x in fold_rows],'n_folds':len(fold_rows),'fold_results':fold_rows,'mean_fold_delta':obs,'sample_sd_fold_delta':float(np.std([x['delta'] for x in fold_rows],ddof=1)) if len(fold_rows)>1 else None,'ci_95':[float(lo),float(hi)],'two_sided_p':float(pv),'contains_zero':bool(lo<=0<=hi),'bootstrap_deltas':arr.tolist()}

def attach_clinical(grid):
 clinical=grid[grid.configuration.eq('Clinical-only')][['sample_key','fold',*MODELS]].rename(columns={m:f'{m}_clinical' for m in MODELS});return grid.merge(clinical,on=['sample_key','fold'],how='left',validate='many_to_one')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--iterations',type=int,default=9999);p.add_argument('--seed',type=int,default=42);p.add_argument('--n-jobs',type=int,default=min(8,os.cpu_count() or 1));a=p.parse_args();repo=Path(__file__).resolve().parents[2];out=repo/'results/phase10';out.mkdir(parents=True,exist_ok=True);grid=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');old=json.loads((repo/'results/phase9/clinical_vs_omics_phase9.json').read_text(encoding='utf-8'))
 tasks=[]
 for cohort in ('mixed','primary'):
  d=attach_clinical(grid[grid.cohort.eq(cohort)].copy())
  for model in MODELS:
   for config in CONFIGS:tasks.append((d,cohort,model,config,a.iterations,a.seed))
 print(f'Running {len(tasks)} comparisons with n_jobs={a.n_jobs}, iterations={a.iterations}, seed={a.seed}',flush=True);rows=Parallel(n_jobs=a.n_jobs,verbose=10)(delayed(one_test)(*x) for x in tasks)
 for cohort in ('mixed','primary'):
  ix=[i for i,r in enumerate(rows) if r['cohort']==cohort];q=multipletests([rows[i]['two_sided_p'] for i in ix],method='fdr_bh')[1]
  for i,v in zip(ix,q):rows[i]['bh_fdr_q']=float(v)
 oldmaps={c:{(r['model'],r['configuration']):r for r in old[c]['comparisons']} for c in ('mixed','primary')};changes=[]
 for r in rows:
  o=oldmaps[r['cohort']][(r['model'],r['configuration'])];r['old_pooled']={'delta':o['delta_c_index'],'ci_95':o['ci_95'],'p':o['bootstrap_p_value'],'bh_fdr_q':o['bh_fdr_p_value'],'contains_zero':o['contains_zero']};r['direction_changed']=bool(np.sign(o['delta_c_index'])!=np.sign(r['mean_fold_delta']));r['raw_ci_significance_changed']=bool(o['contains_zero']!=r['contains_zero']);r['bh_significance_changed']=bool((o['bh_fdr_p_value']<.05)!=(r['bh_fdr_q']<.05))
  if r['direction_changed'] or r['raw_ci_significance_changed'] or r['bh_significance_changed']:changes.append({'cohort':r['cohort'],'model':r['model'],'configuration':r['configuration'],'direction_changed':r['direction_changed'],'raw_ci_significance_changed':r['raw_ci_significance_changed'],'bh_significance_changed':r['bh_significance_changed'],'old':r['old_pooled'],'new':{'delta':r['mean_fold_delta'],'ci_95':r['ci_95'],'p':r['two_sided_p'],'bh_fdr_q':r['bh_fdr_q']}})
 payload={'method':{'estimand':'unweighted mean of estimable outer-fold Harrell-C deltas (configuration minus Clinical-only)','bootstrap':'participant resampling independently within each outer fold','iterations':a.iterations,'seed':a.seed,'bh_family':'24 exploratory comparisons separately within each cohort'},'comparisons':rows,'changed_rows':changes,'deepsurv_retraining_code_audit':{'reinitialized_after_early_stopping':False,'training_data':'internal-training only','preprocessing_refit_after_epoch_selection':False,'outer_test_prediction':'restored best early-stopping checkpoint; no separate retrained model'}};(out/'phase10_exploratory_foldstratified_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
 lines=['# Phase 10 exploratory fold-stratified comparisons and DeepSurv code audit','',f"Method: Harrell C within each estimable outer fold; delta = configuration minus Clinical-only; unweighted mean across estimable folds; participant bootstrap independently within folds; {a.iterations:,} iterations; seed {a.seed}; percentile 95% CI; two-sided empirical tail p-value with +1 finite-sample correction; BH-FDR across 24 tests separately per cohort.",'','## Old pooled versus Phase 10 fold-stratified results','','| Cohort | Model | Configuration | Folds | Old pooled Δ | Old 95% CI | Old p | Old q | New fold-mean Δ | New 95% CI | New p | New q |','|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in rows:
  o=r['old_pooled'];lines.append(f"| {r['cohort']} | {r['model']} | {r['configuration']} | {','.join(map(str,r['folds_used']))} ({r['n_folds']}/5) | {o['delta']:.6f} | [{o['ci_95'][0]:.6f}, {o['ci_95'][1]:.6f}] | {o['p']:.6f} | {o['bh_fdr_q']:.6f} | {r['mean_fold_delta']:.6f} | [{r['ci_95'][0]:.6f}, {r['ci_95'][1]:.6f}] | {r['two_sided_p']:.6f} | {r['bh_fdr_q']:.6f} |")
 lines += ['','## Methylation-only estimable-fold warning','','These rows use each model’s own estimable fold set. Consequently, methylation-only results across CoxPH, RSF, and DeepSurv are not directly comparable.','','| Cohort | Model | Number of folds | Fold IDs |','|---|---|---:|---|']
 for r in rows:
  if r['configuration']=='Methylation only':lines.append(f"| {r['cohort']} | {r['model']} | {r['n_folds']} | {', '.join(map(str,r['folds_used']))} |")
 lines += ['','## Rows with direction or significance changes','']
 if not changes:lines.append('None.')
 for x in changes:lines.append(f"- {x['cohort']} / {x['model']} / {x['configuration']}: direction changed={'yes' if x['direction_changed'] else 'no'}; raw-CI significance changed={'yes' if x['raw_ci_significance_changed'] else 'no'}; BH q<0.05 status changed={'yes' if x['bh_significance_changed'] else 'no'}; old Δ={x['old']['delta']:.6f}, p={x['old']['p']:.6f}, q={x['old']['bh_fdr_q']:.6f}; new Δ={x['new']['delta']:.6f}, p={x['new']['p']:.6f}, q={x['new']['bh_fdr_q']:.6f}.")
 lines += ['','## Code-verified DeepSurv procedure','','The premise that the network is retrained after selecting an epoch is not implemented in the retained Phase 7/9 code. Exact answers:','','1. **No reinitialization/retraining occurs.** `fit_network` trains one network while monitoring internal-validation C-index, stores a copy of the best state dictionary, and after stopping loads that checkpoint (`src/models/phase7_deepsurv.py`, `fit_network`, lines 31–40).','2. **Optimization uses internal-training only.** The loss tensors and feature blocks used by the optimizer are constructed from `xf`/`yf`. Internal validation `xv`/`yv` is used to compute the early-stopping C-index; it is not merged back into training (lines 31–40).','3. **No preprocessing refit follows epoch selection.** Phase 7 fits clinical encoding and methylation imputation on internal train (lines 21–29), and performs internal-train outcome-guided selection before the fit (lines 50–57). Validation/test are transform-only. Phase 9 loads those cached internal-train-derived frames and passes them to the same fit function (`src/models/phase9_deepsurv.py`, lines 31–42).','4. **Outer-test risk comes from the restored best checkpoint.** After `model.load_state_dict(best)`, `predict_risk(model, tb)` produces the outer-test risk. It is not produced by a separately retrained model (`src/models/phase7_deepsurv.py`, line 40).','','Therefore manuscript wording that says the network “is retrained” after epoch selection is not code-accurate. The code performs checkpoint selection and restoration only.','','## Seeds','','- Outer folds: seed 42.','- Internal participant-grouped/event-stratified split: 4200 + outer fold.','- Initial and only DeepSurv training run: 42 + outer fold.','- Phase 10 fold-stratified bootstrap: seed 42, 9,999 iterations.']
 (out/'phase10_exploratory_foldstratified_results.md').write_text('\n'.join(lines),encoding='utf-8')
 csv=[]
 for r in rows:csv.append({'cohort':r['cohort'],'model':r['model'],'configuration':r['configuration'],'folds_used':','.join(map(str,r['folds_used'])),'n_folds':r['n_folds'],'mean_fold_delta':r['mean_fold_delta'],'fold_delta_sample_sd':r['sample_sd_fold_delta'],'ci_95_lower':r['ci_95'][0],'ci_95_upper':r['ci_95'][1],'p_value':r['two_sided_p'],'bh_fdr_q':r['bh_fdr_q'],'contains_zero':r['contains_zero']})
 pd.DataFrame(csv).to_csv(repo/'supplementary/table_s2_exploratory_clinical_vs_omics.csv',index=False);print('Wrote Phase 10 report/JSON and updated Supplementary Table S2')

if __name__=='__main__':main()
