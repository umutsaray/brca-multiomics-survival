#!/usr/bin/env python3
"""Phase 17 Task 7: Uno versus Harrell C-index sensitivity using stored OOF risks."""
from __future__ import annotations
import argparse,json,platform,sys,time
from pathlib import Path
import numpy as np,pandas as pd,sklearn,sksurv
from sksurv.metrics import concordance_index_censored,concordance_index_ipcw
from sksurv.nonparametric import CensoringDistributionEstimator
from sksurv.util import Surv
MODELS=('coxph','rsf','deepsurv');CONFIGS=('Clinical-only','Full');TAUS={'existing_none':None,'3_year':1095.,'5_year':1825.,'7_year':2555.};CONSTRUCTIONS=('other_four_oof','outer_training_outcomes')
def harrell(d,m):return float(concordance_index_censored(d.event.astype(bool),d.time,d[m])[0])
def uno(train,test,m,tau):return float(concordance_index_ipcw(Surv.from_arrays(train.event.astype(bool),train.time),Surv.from_arrays(test.event.astype(bool),test.time),test[m].to_numpy(),tau=tau)[0])
def weights(train,test,tau):
 st=Surv.from_arrays(train.event.astype(bool),train.time);se=Surv.from_arrays(test.event.astype(bool),test.time);base=CensoringDistributionEstimator().fit(st).predict_ipcw(se);mask=np.ones(len(test),bool) if tau is None else test.time.to_numpy()<tau;base=np.where(mask,base,0.);w=base**2;pos=w>0;return {'max_base_ipcw':float(base.max()) if len(base) else None,'max_uno_squared_weight':float(w.max()) if len(w) else None,'effective_sample_size_uno_weights':float(w[pos].sum()**2/(w[pos]@w[pos])) if pos.any() else None,'positive_weight_n':int(pos.sum()),'at_risk_rows_n':None if tau is None else int((test.time>=tau).sum()),'at_risk_participants_n':None if tau is None else int(test.loc[test.time>=tau,'participant'].nunique())}
def safe_weights(train,test,tau):
 try:return {**weights(train,test,tau),'error':None}
 except Exception as e:return {'max_base_ipcw':None,'max_uno_squared_weight':None,'effective_sample_size_uno_weights':None,'positive_weight_n':None,'at_risk_rows_n':None if tau is None else int((test.time>=tau).sum()),'at_risk_participants_n':None if tau is None else int(test.loc[test.time>=tau,'participant'].nunique()),'error':f'{type(e).__name__}: {e}'}
def safe_uno(train,test,m,tau):
 try:return uno(train,test,m,tau),None
 except Exception as e:return None,f'{type(e).__name__}: {e}'
def bootstrap(data,variants,iters,seed):
 rng=np.random.default_rng(seed);arr={(co,model,t,c):np.full(iters,np.nan) for co in data for model in MODELS for t in TAUS for c in CONSTRUCTIONS}
 for co,d in data.items():
  for i in range(iters):
   deltas={(m,t,c):[] for m in MODELS for t in TAUS for c in CONSTRUCTIONS}
   for fold in sorted(d.fold.unique()):
    fd=d[d.fold==fold];base=fd[fd.configuration=='Clinical-only'].sort_values('sample_key').reset_index(drop=True);full=fd[fd.configuration=='Full'].sort_values('sample_key').reset_index(drop=True)
    participants=base.participant.drop_duplicates().to_numpy();chosen=participants[rng.integers(0,len(participants),len(participants))]
    b=pd.concat([base[base.participant==p] for p in chosen],ignore_index=True);f=pd.concat([full[full.participant==p] for p in chosen],ignore_index=True)
    trains={'other_four_oof':d[(d.fold!=fold)&(d.configuration=='Clinical-only')],'outer_training_outcomes':d[(d.fold!=fold)&(d.configuration=='Full')]}
    for m in MODELS:
     for tn,tau in TAUS.items():
      for c,tr in trains.items():
       x,e1=safe_uno(tr,b,m,tau);y,e2=safe_uno(tr,f,m,tau)
       if x is not None and y is not None:deltas[(m,tn,c)].append(x-y)
   for k,v in deltas.items():
    if len(v)==5:arr[(co,*k)][i]=np.mean(v)
  print('bootstrap',co,iters,flush=True)
 out={}
 for (co,m,t,c),v in arr.items():
  good=v[np.isfinite(v)];out.setdefault(co,{}).setdefault(m,{}).setdefault(t,{})[c]={'valid_iterations':len(good),'ci_95':[float(x) for x in np.percentile(good,[2.5,97.5])] if len(good) else [None,None],'bootstrap_deltas':[float(x) if np.isfinite(x) else None for x in v]}
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--iterations',type=int,default=9999);p.add_argument('--seed',type=int,default=42);a=p.parse_args();st=time.perf_counter();repo=a.repo_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);oof=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');oof['participant']=oof.sample_key.astype(str).str[:12];data={c:oof[(oof.cohort==c)&oof.configuration.isin(CONFIGS)].copy() for c in ('mixed','primary')};folds={}
 for co,d in data.items():
  folds[co]={}
  for fold in sorted(d.fold.unique()):
   folds[co][str(fold)]={};test0=d[(d.fold==fold)&(d.configuration=='Clinical-only')].sort_values('sample_key');trains={'other_four_oof':d[(d.fold!=fold)&(d.configuration=='Clinical-only')],'outer_training_outcomes':d[(d.fold!=fold)&(d.configuration=='Full')]};folds[co][str(fold)]['training_constructions_identical_outcomes']=bool(trains['other_four_oof'][['sample_key','event','time']].sort_values('sample_key').reset_index(drop=True).equals(trains['outer_training_outcomes'][['sample_key','event','time']].sort_values('sample_key').reset_index(drop=True)))
   folds[co][str(fold)]['test_participants']=int(test0.participant.nunique());folds[co][str(fold)]['test_rows']=int(len(test0));folds[co][str(fold)]['weight_diagnostics']={c:{tn:safe_weights(tr,test0,tau) for tn,tau in TAUS.items()} for c,tr in trains.items()};folds[co][str(fold)]['metrics']={}
   for cfg in CONFIGS:
    te=d[(d.fold==fold)&(d.configuration==cfg)].sort_values('sample_key');folds[co][str(fold)]['metrics'][cfg]={}
    for m in MODELS:
     rec={'harrell_c':harrell(te,m),'uno':{}}
     for c,tr in trains.items():
      rec['uno'][c]={}
      for tn,tau in TAUS.items():
       value,error=safe_uno(tr,te,m,tau);rec['uno'][c][tn]={'value':value,'error':error}
     folds[co][str(fold)]['metrics'][cfg][m]=rec
 boot=bootstrap(data,None,a.iterations,a.seed);comparisons={}
 for co in data:
  comparisons[co]={}
  for m in MODELS:
   comparisons[co][m]={}
   for tn in TAUS:
    comparisons[co][m][tn]={}
    for c in CONSTRUCTIONS:
     ds=[]
     for f in folds[co].values():
      x=f['metrics']['Clinical-only'][m]['uno'][c][tn]['value'];y=f['metrics']['Full'][m]['uno'][c][tn]['value'];ds.append(None if x is None or y is None else x-y)
     comparisons[co][m][tn][c]={'fold_differences':ds,'mean_difference':float(np.mean(ds)) if all(x is not None for x in ds) else None,**boot[co][m][tn][c]}
 for co in comparisons:
  for m in comparisons[co]:
   for c in CONSTRUCTIONS:
    reference=comparisons[co][m]['existing_none'][c]['mean_difference']
    for tn in TAUS:
     current=comparisons[co][m][tn][c]['mean_difference'];comparisons[co][m][tn][c]['ordering_changed_from_existing_tau']=None if reference is None or current is None else bool(np.sign(reference)!=np.sign(current))
 payload={'task':'Phase 17 Task 7 Uno versus Harrell sensitivity','status':'completed','decision_rules_fixed_before_results':{'stored_predictions_only':True,'participant_bootstrap_key':'first 12 characters of sample_key','cohorts':['mixed','primary'],'configurations':list(CONFIGS),'models':list(MODELS),'taus_days':TAUS,'weight_constructions':list(CONSTRUCTIONS),'delta':'Clinical-only minus Full','bootstrap_iterations':a.iterations,'bootstrap_seed':a.seed},'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__,'sksurv':sksurv.__version__},'command':' '.join([sys.executable,*sys.argv]),'runtime_seconds':time.perf_counter()-st,'folds':folds,'comparisons':comparisons};(out/'task7_uno_harrell_sensitivity.json').write_text(json.dumps(payload,indent=2),encoding='utf8');lines=['# Phase 17 Task 7 — Uno versus Harrell C-index sensitivity','',f"Design fixed before results: `{json.dumps(payload['decision_rules_fixed_before_results'])}`",'','## Per-fold Harrell and Uno C-index','','| Cohort | Fold | Configuration | Model | Harrell | Uno existing tau | Uno 3-year | Uno 5-year | Uno 7-year |','|---|---:|---|---|---:|---:|---:|---:|---:|']
 for co,fs in folds.items():
  for fold,f in fs.items():
   for cfg,models in f['metrics'].items():
    for m,z in models.items():
     u=z['uno']['other_four_oof'];lines.append(f"| {co} | {fold} | {cfg} | {m} | {z['harrell_c']} | {u['existing_none']['value']} | {u['3_year']['value']} | {u['5_year']['value']} | {u['7_year']['value']} |")
 lines+=['','## IPCW weight diagnostics and participants at risk','','| Cohort | Fold | Tau | Weight construction | Participants at risk | Maximum IPCW | Maximum squared Uno weight | Weight ESS |','|---|---:|---|---|---:|---:|---:|---:|']
 for co,fs in folds.items():
  for fold,f in fs.items():
   for c,ts in f['weight_diagnostics'].items():
    for tn,z in ts.items():lines.append(f"| {co} | {fold} | {tn} | {c} | {z['at_risk_participants_n']} | {z['max_base_ipcw']} | {z['max_uno_squared_weight']} | {z['effective_sample_size_uno_weights']} |")
 lines+=['','## Clinical-only minus Full Uno differences','','| Cohort | Model | Tau | Weight construction | Mean difference | 95% CI | Valid bootstrap iterations | Ordering changed vs existing tau |','|---|---|---|---|---:|---|---:|---|']
 for co,ms in comparisons.items():
  for m,ts in ms.items():
   for tn,cs in ts.items():
    for c,z in cs.items():lines.append(f"| {co} | {m} | {tn} | {c} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['valid_iterations']} | {z['ordering_changed_from_existing_tau']} |")
 lines+=['','## Weight diagnostics','','Maximum base IPCW, maximum squared Uno weight, squared-weight effective sample size, positive-weight count, at-risk count, and any error are retained per cohort/fold/tau/construction without rounding in the JSON. The `other_four_oof` path uses outcomes from the other four OOF folds; the independently assembled `outer_training_outcomes` path uses the corresponding Full-configuration outer-training outcome rows. Their equality is checked and recorded for every fold.','',f"Runtime seconds: `{payload['runtime_seconds']}`."];(out/'task7_uno_harrell_sensitivity.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
if __name__=='__main__':main()
