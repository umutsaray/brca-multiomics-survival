#!/usr/bin/env python3
"""Finalize/render a completed Phase 17 Task 7 JSON without recomputing statistics."""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np,pandas as pd

TAUS={'existing_none':None,'3_year':1095.,'5_year':1825.,'7_year':2555.}
CONSTRUCTIONS=('other_four_oof','outer_training_outcomes')

def finite_json(x):
 if isinstance(x,float) and not math.isfinite(x):return None
 if isinstance(x,list):return [finite_json(v) for v in x]
 if isinstance(x,dict):return {k:finite_json(v) for k,v in x.items()}
 return x

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-dir',type=Path,required=True);p.add_argument('--json',type=Path,required=True);p.add_argument('--markdown',type=Path,required=True);a=p.parse_args();repo=a.repo_dir.resolve();jp=a.json.resolve();mp=a.markdown.resolve();payload=json.loads(jp.read_text(encoding='utf-8'));oof=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');oof['participant']=oof.sample_key.astype(str).str[:12]
 for co,fs in payload['folds'].items():
  for fold,f in fs.items():
   test=oof[(oof.cohort==co)&(oof.fold==int(fold))&(oof.configuration=='Clinical-only')]
   for construction,ts in f['weight_diagnostics'].items():
    for tn,tau in TAUS.items():
     z=ts[tn];z.pop('at_risk_n',None);z['at_risk_rows_n']=None if tau is None else int((test.time>=tau).sum());z['at_risk_participants_n']=None if tau is None else int(test.loc[test.time>=tau,'participant'].nunique())
 for co in payload['comparisons']:
  for model in payload['comparisons'][co]:
   for construction in CONSTRUCTIONS:
    reference=payload['comparisons'][co][model]['existing_none'][construction]['mean_difference']
    for tn in TAUS:
     z=payload['comparisons'][co][model][tn][construction];current=z['mean_difference'];z.pop('ordering_changed_from_positive',None);z['ordering_changed_from_existing_tau']=None if reference is None or current is None else bool(np.sign(reference)!=np.sign(current))
 payload['report_finalization']={'command':' '.join([sys.executable,*sys.argv]),'statistics_recomputed':False,'derived_fields':['at-risk participant counts from stored OOF outcomes','ordering-change booleans from stored mean differences','non-finite bootstrap draws serialized as null'],'reason':'Make participant counting explicit and JSON standards-compliant; no C-index, confidence interval, or weight statistic was recomputed.'}
 payload=finite_json(payload);jp.write_text(json.dumps(payload,indent=2,allow_nan=False),encoding='utf-8')
 lines=['# Phase 17 Task 7 — Uno versus Harrell C-index sensitivity','','**Exploratory sensitivity analysis.** No model was refit; only stored out-of-fold predictions were used.','',f"Design fixed before results: `{json.dumps(payload['decision_rules_fixed_before_results'])}`",'',f"Analysis runtime seconds: `{payload['runtime_seconds']}`.",'','## Per-fold Harrell and Uno C-index','','The Uno columns below use the existing censoring-weight construction (outcomes from the other four OOF folds). Values for the independently assembled outer-training-outcome construction are retained in the JSON.','','| Cohort | Fold | Configuration | Model | Harrell | Uno existing tau | Uno 3-year | Uno 5-year | Uno 7-year |','|---|---:|---|---|---:|---:|---:|---:|---:|']
 for co,fs in payload['folds'].items():
  for fold,f in fs.items():
   for cfg,models in f['metrics'].items():
    for model,z in models.items():
     u=z['uno']['other_four_oof'];lines.append(f"| {co} | {fold} | {cfg} | {model} | {z['harrell_c']} | {u['existing_none']['value']} | {u['3_year']['value']} | {u['5_year']['value']} | {u['7_year']['value']} |")
 lines+=['','## IPCW diagnostics and participants at risk','','`Maximum IPCW` is the maximum base event weight 1/G; `maximum squared Uno weight` is the corresponding weight used by Uno concordance. ESS is computed from the squared Uno weights. Participant counts use the 12-character TCGA barcode.','','| Cohort | Fold | Tau | Weight construction | Participants at risk | Rows at risk | Maximum IPCW | Maximum squared Uno weight | Weight ESS | Error |','|---|---:|---|---|---:|---:|---:|---:|---:|---|']
 for co,fs in payload['folds'].items():
  for fold,f in fs.items():
   for construction,ts in f['weight_diagnostics'].items():
    for tn,z in ts.items():lines.append(f"| {co} | {fold} | {tn} | {construction} | {z['at_risk_participants_n']} | {z['at_risk_rows_n']} | {z['max_base_ipcw']} | {z['max_uno_squared_weight']} | {z['effective_sample_size_uno_weights']} | {z['error']} |")
 lines+=['','## Clinical-only minus Full Uno differences','','The estimand is the unweighted mean of five within-fold C-index differences. Positive values favor Clinical-only. The bootstrap resamples 12-character participants independently within each fold. A valid-iteration count below 9,999 means that at least one resampled fold had no estimable comparable pair; such draws are stored as `null`, not silently replaced. Ordering change compares each tau with the existing-tau sign under the same weight construction.','','| Cohort | Model | Tau | Weight construction | Mean difference | 95% CI | Valid / requested | Ordering changed vs existing tau |','|---|---|---|---|---:|---|---:|---|']
 for co,models in payload['comparisons'].items():
  for model,taus in models.items():
   for tn,constructions in taus.items():
    for construction,z in constructions.items():lines.append(f"| {co} | {model} | {tn} | {construction} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['valid_iterations']} / {payload['decision_rules_fixed_before_results']['bootstrap_iterations']} | {z['ordering_changed_from_existing_tau']} |")
 lines+=['','## Construction audit','']
 for co,fs in payload['folds'].items():
  lines.append(f"- {co}: outer-training outcomes and other-four-OOF outcomes were identical in all folds = `{all(f['training_constructions_identical_outcomes'] for f in fs.values())}`.")
 errors=[]
 for co,fs in payload['folds'].items():
  for fold,f in fs.items():
   for construction,ts in f['weight_diagnostics'].items():
    for tn,z in ts.items():
     if z['error'] is not None:errors.append([co,fold,construction,tn,z['error']])
 lines+=['',f"Fold-level IPCW diagnostic errors: `{json.dumps(errors)}`.",'','All raw bootstrap draws, fold-level values, row/participant counts, and weight diagnostics are retained without rounding in the JSON.','',f"Report finalization: `{json.dumps(payload['report_finalization'])}`."]
 mp.write_text('\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':main()
