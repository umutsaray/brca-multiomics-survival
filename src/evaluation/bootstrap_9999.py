#!/usr/bin/env python3
"""Repeat the pre-specified Full-minus-Clinical paired bootstrap 9,999 times."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored

MODELS=('coxph','rsf','deepsurv')
def ci(d,c):return float(concordance_index_censored(d.event.astype(bool),d.time,d[c])[0])
def full_test(data,iterations,seed):
 w=data[data.configuration.isin(['Clinical-only','Full'])].pivot(index=['sample_key','fold','event','time'],columns='configuration',values=list(MODELS));w.columns=[f'{m}::{c}' for m,c in w.columns];w=w.reset_index();pids=w.sample_key.str[:12].to_numpy();unique=pd.unique(pids);rng=np.random.default_rng(seed);samples=[]
 for _ in range(iterations):
  chosen=rng.choice(unique,size=len(unique),replace=True);samples.append(np.concatenate([np.flatnonzero(pids==p) for p in chosen]))
 results=[];raw={}
 for m in MODELS:
  diffs=np.asarray([ci(w.iloc[idx],f'{m}::Full')-ci(w.iloc[idx],f'{m}::Clinical-only') for idx in samples if w.iloc[idx].event.sum()>0]);lo,hi=np.percentile(diffs,[2.5,97.5]);p=min(1.,2*min((np.count_nonzero(diffs<=0)+1)/(len(diffs)+1),(np.count_nonzero(diffs>=0)+1)/(len(diffs)+1)));point=ci(w,f'{m}::Full')-ci(w,f'{m}::Clinical-only');results.append({'model':m,'delta_c_index':point,'ci_95':[float(lo),float(hi)],'bootstrap_p_value':float(p),'valid_bootstrap_iterations':len(diffs)});raw[m]=diffs.tolist()
 return {'n_rows':len(w),'n_participants':len(unique),'iterations_requested':iterations,'seed':seed,'comparisons':results,'raw_paired_differences':raw}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--iterations',type=int,default=9999);p.add_argument('--random-state',type=int,default=42);p.add_argument('--output',type=Path,default=Path('bootstrap_9999_primary.md'));a=p.parse_args();payload={}
 for cohort,file in (('mixed','ablation_oof_extended.parquet'),('primary','ablation_oof_primary.parquet')):
  payload[cohort]=full_test(pd.read_parquet(a.input_dir/file),a.iterations,a.random_state)
 lines=['# 9,999-iteration paired participant bootstrap','',f'- Iterations requested: {a.iterations}',f'- Random seed: {a.random_state}','- Resampling unit: 12-character TCGA participant ID','- Comparison direction: Full multi-omics minus Clinical-only','- Input risks: unchanged stored out-of-fold predictions','']
 for cohort,x in payload.items():
  lines += [f'## {cohort}','','| Model | Δ C-index | 95% CI | Two-sided bootstrap p-value | Valid iterations |','|---|---:|---|---:|---:|']
  for r in x['comparisons']:lines.append(f"| {r['model']} | {r['delta_c_index']} | [{r['ci_95'][0]}, {r['ci_95'][1]}] | {r['bootstrap_p_value']} | {r['valid_bootstrap_iterations']} |")
  lines.append('')
 a.output.write_text('\n'.join(lines),encoding='utf-8');a.output.with_suffix('.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');print(f'Wrote {a.output} and {a.output.with_suffix(".json")}')
if __name__=='__main__':main()
