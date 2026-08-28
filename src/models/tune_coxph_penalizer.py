#!/usr/bin/env python3
"""Nested-CV tuning of the CoxPH ridge penalizer."""
import argparse,json,time,warnings
from pathlib import Path
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sksurv.metrics import concordance_index_censored
from phase1_common import load_cohort,grouped_folds,prepare_fold,rank_fold,select_from_rankings

def fit_score(xtr,ytr,xte,yte,p,max_steps=100):
 ev=ytr.event.astype(bool);a=xtr.loc[ev].var();b=xtr.loc[~ev].var();cols=[c for c in xtr if a[c]>1e-12 and b[c]>1e-12];sc=StandardScaler();u=sc.fit_transform(xtr[cols]);v=sc.transform(xte[cols]);d=pd.DataFrame(u,columns=cols);d['time']=ytr.time.to_numpy();d['event']=ytr.event.to_numpy()
 with warnings.catch_warnings(record=True) as w:
  warnings.simplefilter('always');m=CoxPHFitter(penalizer=p);m.fit(d,'time','event',fit_options={'max_steps':max_steps})
 risk=m.predict_partial_hazard(pd.DataFrame(v,columns=cols)).to_numpy().ravel();return float(concordance_index_censored(yte.event.astype(bool),yte.time,risk)[0]),[str(z.message) for z in w]
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--cohort',choices=['mixed','primary'],default='mixed');p.add_argument('--penalizers',type=float,nargs='+',default=[.01,.05,.1,.5,1.]);p.add_argument('--top-n',type=int,default=200);p.add_argument('--random-state',type=int,default=42);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--methylation-chunksize',type=int,default=5000);a=p.parse_args();root=a.input_dir.resolve();base,y=load_cohort(root,a.cohort);suffix='_primary' if a.cohort=='primary' else '';outer=[];conv={str(v):0 for v in a.penalizers};start=time.perf_counter()
 for fold,(tr,te) in enumerate(grouped_folds(y,5,a.random_state),1):
  xtr,xte,_=prepare_fold(root,base,y,tr,te,a.meth_prefilter,a.knn_neighbors,a.methylation_chunksize);ytr,yte=y.iloc[tr],y.iloc[te];rank=rank_fold(xtr,ytr);sel,_=select_from_rankings(rank,xtr,a.top_n,True);xtr=xtr[sel];xte=xte[sel];groups=pd.Index(ytr.index).str[:12];inner=StratifiedGroupKFold(3,shuffle=True,random_state=a.random_state+fold);scores={v:[] for v in a.penalizers}
  for it,iv in inner.split(np.zeros(len(ytr)),ytr.event,groups):
   for v in a.penalizers:
    s,w=fit_score(xtr.iloc[it],ytr.iloc[it],xtr.iloc[iv],ytr.iloc[iv],v);scores[v].append(s);conv[str(v)]+=int(bool(w))
  means={v:float(np.mean(z)) for v,z in scores.items()};best=max(a.penalizers,key=lambda v:(means[v],-v));test,w=fit_score(xtr,ytr,xte,yte,best);conv[str(best)]+=int(bool(w));outer.append({'fold':fold,'selected_penalizer':best,'inner_means':means,'test_cindex':test,'outer_warning':bool(w)});print(f'Fold {fold}: selected={best!r}, test_cindex={test!r}, warnings={w}')
 vals=[z['test_cindex'] for z in outer];payload={'cohort':a.cohort,'folds':outer,'mean':float(np.mean(vals)),'sample_sd':float(np.std(vals,ddof=1)),'convergence_warning_fit_counts':conv,'elapsed_seconds':time.perf_counter()-start};(root/f'coxph_penalizer_tuning{suffix}.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');(root/f'coxph_penalizer_tuning{suffix}.md').write_text(f"# Nested CoxPH Penalizer Tuning\n\nOuter C-index: {payload['mean']} ± {payload['sample_sd']}\n\nConvergence-warning fit counts: `{conv}`\n",encoding='utf-8')
if __name__=='__main__':main()
