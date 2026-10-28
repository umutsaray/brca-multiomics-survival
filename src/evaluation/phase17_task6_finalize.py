#!/usr/bin/env python3
"""Finalize a stopped Phase 17 Task 6 report without fitting or resampling."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path

CONFIGS=('Clinical-only','Full')

def main():
 p=argparse.ArgumentParser();p.add_argument('--json',type=Path,required=True);p.add_argument('--markdown',type=Path,required=True);p.add_argument('--checkpoint-dir',type=Path,required=True);a=p.parse_args();jp=a.json.resolve();mp=a.markdown.resolve();payload=json.loads(jp.read_text(encoding='utf-8'));payload['completed_checkpoint_count']=len(list(a.checkpoint_dir.resolve().glob('fold*_fraction*_rep*.json')))
 for fraction,z in payload['clinical_minus_full_bootstrap'].items():
  if z is not None:
   s=payload['summary'][fraction];z['point_difference']=s['Clinical-only']['mean']-s['Full']['mean']
 payload['partial_result_warning']='The 50% values contain one completed outer fold and one repetition only, not the pre-specified five-fold repetition-averaged estimand; the 75% analysis was not started. Neither is interpreted as a completed learning-curve result.'
 payload['report_finalization']={'command':' '.join([sys.executable,*sys.argv]),'model_fits_run':0,'bootstrap_resamples_run':0,'derived_fields':['correct checkpoint-file count','stored-summary Clinical-only minus Full point difference'],'statistics_replaced':False}
 jp.write_text(json.dumps(payload,indent=2,allow_nan=False),encoding='utf-8')
 S=payload['summary'];B=payload['clinical_minus_full_bootstrap'];lines=['# Phase 17 Task 6 — Reduced CoxPH-only learning curve','','Status: **stopped**.','',payload['partial_result_warning'],'','RSF and DeepSurv were skipped as a pre-fixed runtime decision. CoxPH used the Table 7 fixed ridge penalizer 0.1. Fractions 50% and 75% were specified for three repetitions; 100% is the stored Table 7 result.','',f"Stop reason: {payload['stop_reason']}",f"Projected remaining seconds at stop: `{payload['projected_remaining_seconds_at_stop']}`.",f"Completed checkpoint count: `{payload['completed_checkpoint_count']}`.",'','## C-index summary','','| Fraction | Configuration | Completed folds | Mean | Sample SD | Fold values |','|---:|---|---:|---:|---:|---|']
 for fraction in ('0.5','0.75','1.0'):
  for config in CONFIGS:
   z=S[fraction][config];lines.append(f"| {fraction} | {config} | {z['completed_folds']} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_repetition_mean_cindices'])}` |")
 lines+=['','## Completed reduced-fit checkpoint','','| Fraction | Fold | Repetition | Configuration | Training participants | Training events | Selected RNA | Selected miRNA | Selected methylation | Selected clinical | C-index |','|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|']
 for config in CONFIGS:
  z=S['0.5'][config]['folds'][0];c=z['selected_features_per_block'][0];lines.append(f"| 0.5 | {z['fold']} | 1 | {config} | {z['training_participants'][0]} | {z['training_events'][0]} | {c['rna']} | {c['mirna']} | {c['methylation']} | {c['clinical']} | {z['repetition_cindices'][0]} |")
 lines+=['','## Clinical-only minus Full fold-stratified paired bootstrap','','| Fraction | Point difference | 95% CI | Status |','|---:|---:|---|---|']
 for fraction in ('0.5','0.75','1.0'):
  z=B[fraction];lines.append(f"| {fraction} | {None if z is None else z['point_difference']} | `{json.dumps(None if z is None else z['ci_95'])}` | {'unavailable: incomplete design' if z is None else 'available'} |")
 lines+=['','## Observed trend','','No 50%–75%–100% trend classification is reported because the 50% design is incomplete and the 75% design was not started. No extrapolation was performed.','',f"Actual elapsed time (seconds): `{payload['runtime_seconds']}`.",'',f"Report finalization: `{json.dumps(payload['report_finalization'])}`."]
 mp.write_text('\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':main()
