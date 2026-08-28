#!/usr/bin/env python3
"""Patient-level paired bootstrap comparisons on pooled out-of-fold risks.

The statistic is the difference in Harrell C-index on paired OOF predictions.
Each iteration samples 12-character TCGA participant IDs with replacement and
includes all records belonging to every sampled ID.  The same sampled row
indices are used for all three models.  Percentile limits use the 2.5th and
97.5th percentiles.  The two-sided bootstrap tail p-value is twice the smaller
finite-sample-corrected proportion of differences <=0 or >=0.
"""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from sksurv.metrics import concordance_index_censored
MODELS=['coxph','rsf','deepsurv']
def ci(e,t,r):return float(concordance_index_censored(e.astype(bool),t,r)[0])
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--cohort',choices=['mixed','primary'],default='mixed');p.add_argument('--iterations',type=int,default=1000);p.add_argument('--random-state',type=int,default=42);a=p.parse_args();root=a.input_dir.resolve();suffix='_primary' if a.cohort=='primary' else '';d=pd.read_parquet(root/f'phase1_oof_predictions{suffix}.parquet');groups=d.sample_key.str[:12];patients=groups.unique();rng=np.random.default_rng(a.random_state);raw={m:[] for m in MODELS}
 for _ in range(a.iterations):
  chosen=rng.choice(patients,len(patients),replace=True);idx=np.concatenate([np.flatnonzero(groups.to_numpy()==g) for g in chosen]);
  if d.event.iloc[idx].sum()==0:continue
  for m in MODELS:raw[m].append(ci(d.event.to_numpy()[idx],d.time.to_numpy()[idx],d[m].to_numpy()[idx]))
 point={m:ci(d.event,d.time,d[m]) for m in MODELS};comparisons={}
 for other in ['coxph','rsf']:
  diff=np.asarray(raw['deepsurv'])-np.asarray(raw[other]);lo,hi=np.percentile(diff,[2.5,97.5]);pval=float(2*min((np.sum(diff<=0)+1)/(len(diff)+1),(np.sum(diff>=0)+1)/(len(diff)+1)));comparisons[f'deepsurv_vs_{other}']={'point_difference':point['deepsurv']-point[other],'ci_95':[float(lo),float(hi)],'contains_zero':bool(lo<=0<=hi),'paired_bootstrap_p_value':min(1.,pval),'interpretation':'anlamlı fark bulunamadı' if lo<=0<=hi else 'bootstrap CI sıfırı içermiyor'}
 payload={'cohort':a.cohort,'n_rows':len(d),'n_participants':len(patients),'iterations_requested':a.iterations,'iterations_valid':len(raw['coxph']),'oof_cindex':point,'comparisons':comparisons,'raw_bootstrap_cindices':raw};(root/f'statistical_comparison_v2{suffix}.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');lines=['# OOF Paired Statistical Comparison','',f"Cohort: {a.cohort}; rows: {len(d)}; participants: {len(patients)}; valid bootstrap: {len(raw['coxph'])}.",'', '| Comparison | Difference | 95% CI | p-value | Interpretation |','|---|---:|---:|---:|---|'];
 for k,v in comparisons.items():lines.append(f"| {k} | {v['point_difference']} | [{v['ci_95'][0]}, {v['ci_95'][1]}] | {v['paired_bootstrap_p_value']} | {v['interpretation']} |")
 (root/f'statistical_comparison_v2{suffix}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print(json.dumps(payload|{'raw_bootstrap_cindices':'written to JSON'},indent=2))
if __name__=='__main__':main()
