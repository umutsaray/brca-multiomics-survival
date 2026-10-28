"""Shared, leakage-safe helpers for Phase 16 exploratory sensitivity analyses."""
from __future__ import annotations
import warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from statsmodels.stats.multitest import multipletests

RIDGE_GRID=(0.01,0.05,0.1,0.5,1.0)

def ci(y,r): return float(concordance_index_censored(y.event.astype(bool),y.time,np.asarray(r,float))[0])
def inner_folds(y,seed,n=5):
 g=y.index.str[:12];return list(StratifiedGroupKFold(n_splits=n,shuffle=True,random_state=seed).split(np.zeros(len(y)),y.event,g))
def diagnostics(a,w):
 rank=int(np.linalg.matrix_rank(a));return {'n_rows':int(a.shape[0]),'n_columns':int(a.shape[1]),'matrix_rank':rank,'rank_deficient':rank<a.shape[1],'convergence_warnings':[str(x.message) for x in w]}
def cox_fit_predict(xtr,ytr,xte,p,max_steps=100):
 ev=ytr.event.astype(bool);v1=xtr.loc[ev].var();v0=xtr.loc[~ev].var();cols=[c for c in xtr if v1[c]>1e-12 and v0[c]>1e-12]
 sc=StandardScaler();a=sc.fit_transform(xtr[cols]);b=sc.transform(xte[cols]);d=pd.DataFrame(a,columns=cols);d['time']=ytr.time.to_numpy();d['event']=ytr.event.to_numpy()
 with warnings.catch_warnings(record=True) as w:
  warnings.simplefilter('always');m=CoxPHFitter(penalizer=p);m.fit(d,'time','event',fit_options={'max_steps':max_steps})
 r=m.predict_log_partial_hazard(pd.DataFrame(b,columns=cols)).to_numpy().ravel();return r,diagnostics(a,w)
def tune_cox(x,y,folds,grid=RIDGE_GRID,max_steps=100):
 scores={str(p):[] for p in grid}
 for tr,va in folds:
  for p in grid:
   r,_=cox_fit_predict(x.iloc[tr],y.iloc[tr],x.iloc[va],p,max_steps);scores[str(p)].append(ci(y.iloc[va],r))
 means={k:float(np.mean(v)) for k,v in scores.items()};best=max(grid,key=lambda p:(means[str(p)],-p));return float(best),scores,means
def rsf_predict(xtr,ytr,xte,seed,trees=500):
 m=RandomSurvivalForest(n_estimators=trees,min_samples_split=10,min_samples_leaf=5,random_state=seed,n_jobs=-1);m.fit(xtr.to_numpy('float32'),Surv.from_arrays(ytr.event.astype(bool),ytr.time));return m.predict(xte.to_numpy('float32'))
def summarize(scores):
 return {k:{'fold_cindices':v,'mean':float(np.mean(v)),'sample_sd':float(np.std(v,ddof=1))} for k,v in scores.items()}
def paired_bootstrap(fold_frames,comparisons,iterations=9999,seed=42):
 rng=np.random.default_rng(seed);arrays={(name,model):np.empty(iterations) for name,left,right in comparisons for model in set(fold_frames[0]['risks'][left])&set(fold_frames[0]['risks'][right])}
 for i in range(iterations):
  delta={k:[] for k in arrays}
  for f in fold_frames:
   y=f['y'];n=len(y);draw=rng.integers(0,n,n)
   vals={(cfg,m):ci(y.iloc[draw],np.asarray(r)[draw]) for cfg,models in f['risks'].items() for m,r in models.items()}
   for name,left,right in comparisons:
    for m in set(f['risks'][left])&set(f['risks'][right]):delta[(name,m)].append(vals[(left,m)]-vals[(right,m)])
  for k in arrays:arrays[k][i]=np.mean(delta[k])
 out={}
 for name,left,right in comparisons:
  out[name]={}
  for m in set(fold_frames[0]['risks'][left])&set(fold_frames[0]['risks'][right]):
   obs=[ci(f['y'],f['risks'][left][m])-ci(f['y'],f['risks'][right][m]) for f in fold_frames];a=arrays[(name,m)];lo,hi=np.percentile(a,[2.5,97.5]);pl=(np.count_nonzero(a<=0)+1)/(iterations+1);ph=(np.count_nonzero(a>=0)+1)/(iterations+1)
   out[name][m]={'fold_differences':obs,'mean_difference':float(np.mean(obs)),'ci_95':[float(lo),float(hi)],'raw_two_sided_p':float(min(1,2*min(pl,ph))),'bootstrap_differences':a.tolist()}
 keys=[(n,m) for n,_,_ in comparisons for m in out[n]];q=multipletests([out[n][m]['raw_two_sided_p'] for n,m in keys],method='fdr_bh')[1]
 for (n,m),v in zip(keys,q):out[n][m]['bh_fdr_q_within_family']=float(v)
 return out
def jaccard_sets(sets):
 vals=[]
 for i in range(len(sets)):
  for j in range(i+1,len(sets)):
   a,b=set(sets[i]),set(sets[j]);vals.append({'folds':[i+1,j+1],'intersection':len(a&b),'union':len(a|b),'jaccard':len(a&b)/len(a|b) if a|b else None})
 return {'pairs':vals,'mean_pairwise_jaccard':float(np.mean([x['jaccard'] for x in vals if x['jaccard'] is not None])) if any(x['jaccard'] is not None for x in vals) else None}
