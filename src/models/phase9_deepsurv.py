#!/usr/bin/env python3
"""Retrain the remaining ablation DeepSurv configurations using Phase 7 preprocessing."""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch

HERE=Path(__file__).resolve().parent
for sub in ("common","data_prep","preprocessing","feature_selection","models","evaluation","explainability"):
    sys.path.insert(0,str(HERE.parents[0]/sub))
from phase1_common import load_cohort, grouped_folds  # noqa: E402
from phase7_deepsurv import inner_split, fit_network  # noqa: E402

COMBOS={
 "RNA-seq only":["rna_"],"Methylation only":["meth_"],"miRNA only":["mirna_"],
 "Clinical + RNA-seq":["clin_","rna_"],
 "Clinical + RNA-seq + Methylation":["clin_","rna_","meth_"],
 "Clinical + RNA-seq + miRNA":["clin_","rna_","mirna_"],
 "Clinical + Methylation":["clin_","meth_"],"Clinical + miRNA":["clin_","mirna_"],
}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--random-state',type=int,default=42);p.add_argument('--epochs',type=int,default=300);p.add_argument('--patience',type=int,default=30);p.add_argument('--learning-rate',type=float,default=.001);p.add_argument('--weight-decay',type=float,default=.0001);p.add_argument('--dropout',type=float,default=.3);a=p.parse_args()
 torch.set_num_threads(1);root=a.input_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);started=time.perf_counter();oof=[]
 payload={'method':{'source':'Phase 7 internal-train-only cached preprocessing/selection','outer_seed':a.random_state,'internal_seed_schema':'4200 + outer_fold','training_seed_schema':f'{a.random_state} + outer_fold','epochs':a.epochs,'patience':a.patience,'learning_rate':a.learning_rate,'weight_decay':a.weight_decay,'dropout':a.dropout,'configurations':list(COMBOS)},'cohorts':{}}
 for cohort in ('mixed','primary'):
  _,y=load_cohort(root,cohort);fold_results=[]
  for fold,(otr,ote) in enumerate(grouped_folds(y,5,a.random_state),1):
   iseed=4200+fold;tseed=a.random_state+fold;fit_idx,val_idx=inner_split(y,otr,iseed);yf,yv,yt=y.iloc[fit_idx],y.iloc[val_idx],y.iloc[ote]
   stem=f'{cohort}_fold{fold}_innerseed{iseed}';cache=root/'phase7_inner_cache';fp,vp,tp,mp=[cache/f'{stem}_{x}' for x in ('fit.parquet','validation.parquet','test.parquet','meta.json')]
   missing=[str(x) for x in (fp,vp,tp,mp) if not x.is_file()]
   if missing: raise FileNotFoundError('Phase 7 cache is required; missing: '+', '.join(missing))
   xf=pd.read_parquet(fp);xv=pd.read_parquet(vp);xt=pd.read_parquet(tp);meta=json.loads(mp.read_text(encoding='utf-8'))
   fr={'fold':fold,'outer_train_n':len(otr),'internal_train_n':len(fit_idx),'internal_validation_n':len(val_idx),'outer_test_n':len(ote),'outer_seed':a.random_state,'internal_split_seed':iseed,'training_seed':tseed,'selection':meta['selection'],'configurations':{}}
   for config,prefixes in COMBOS.items():
    cols=[c for c in xf.columns if c.startswith(tuple(prefixes))]
    if not cols:
     fr['configurations'][config]={'status':'skipped','reason':'zero FDR-qualified features for the requested single block','n_features':0}
     print(cohort,fold,config,'SKIPPED: zero features',flush=True);continue
    result=fit_network(xf[cols],yf,xv[cols],yv,xt[cols],yt,tseed,a.epochs,a.patience,a.learning_rate,a.weight_decay,a.dropout)
    fr['configurations'][config]={'status':'fitted','n_features':len(cols),'cindex':result['cindex'],'best_epoch':result['best_epoch'],'validation_cindex':result['validation_cindex']}
    print(cohort,fold,config,fr['configurations'][config],flush=True)
    for j,key in enumerate(yt.index):oof.append({'sample_key':key,'cohort':cohort,'fold':fold,'configuration':config,'event':int(yt.iloc[j].event),'time':float(yt.iloc[j].time),'deepsurv':float(result['risk'][j])})
   fold_results.append(fr)
  payload['cohorts'][cohort]={'n_rows':len(y),'n_participants':int(y.index.str[:12].nunique()),'folds':fold_results}
 payload['fits_requested']=2*len(COMBOS)*5;payload['fits_completed']=sum(x['status']=='fitted' for c in payload['cohorts'].values() for f in c['folds'] for x in f['configurations'].values());payload['fits_skipped']=payload['fits_requested']-payload['fits_completed'];payload['elapsed_seconds']=time.perf_counter()-started
 pd.DataFrame(oof).to_parquet(out/'phase9_deepsurv_oof.parquet',index=False);(out/'phase9_deepsurv_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');print(f"Completed {payload['fits_completed']}/{payload['fits_requested']} fits; skipped {payload['fits_skipped']}; wrote {out}")

if __name__=='__main__':main()
