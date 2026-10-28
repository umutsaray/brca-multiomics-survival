#!/usr/bin/env python3
"""Phase 16 Task 4: alternative feature-selection sensitivity analyses."""
from __future__ import annotations
import argparse,json,platform,sys,time,warnings
from pathlib import Path
import lifelines,numpy as np,pandas as pd,sklearn,sksurv
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
HERE=Path(__file__).resolve().parent
for sub in ('common','data_prep','preprocessing','feature_selection','models','evaluation'):sys.path.insert(0,str(HERE.parents[0]/sub))
from phase1_common import load_cohort,grouped_folds,prepare_fold
from phase16_common import RIDGE_GRID,ci,cox_fit_predict,inner_folds,tune_cox,rsf_predict,summarize,paired_bootstrap,jaccard_sets
PREFIXES={'rna':'rna_','mirna':'mirna_','methylation':'meth_'};ALPHAS=(0.0001,0.001,0.01,0.1,1.0)

def variance_select(x,n):
 out=[];rank={}
 for b,p in PREFIXES.items():
  s=x[[c for c in x if c.startswith(p)]].var().sort_values(ascending=False);rank[b]=s.index.tolist();out+=rank[b][:min(n,len(s))]
 out += [c for c in x if c.startswith('clin_')];return out,rank
def coxnet_fit(xtr,ytr,xte,alpha):
 sc=StandardScaler();a=sc.fit_transform(xtr);b=sc.transform(xte);m=CoxnetSurvivalAnalysis(l1_ratio=.5,alphas=[alpha],max_iter=100000,tol=1e-7);m.fit(a,Surv.from_arrays(ytr.event.astype(bool),ytr.time));return m.predict(b),m.coef_[:,0],sc
def tune_enet(x,y,folds):
 scores={str(v):[] for v in ALPHAS}
 for tr,va in folds:
  sc=StandardScaler();u=sc.fit_transform(x.iloc[tr]);v=sc.transform(x.iloc[va])
  try:
   model=CoxnetSurvivalAnalysis(l1_ratio=.5,alphas=list(reversed(ALPHAS)),max_iter=100000,tol=1e-7);model.fit(u,Surv.from_arrays(y.iloc[tr].event.astype(bool),y.iloc[tr].time))
   for alpha in ALPHAS:scores[str(alpha)].append(ci(y.iloc[va],model.predict(v,alpha=alpha)))
  except Exception:
   for alpha in ALPHAS:scores[str(alpha)].append(float('-inf'))
 means={k:float(np.mean(z)) for k,z in scores.items()};best=max(ALPHAS,key=lambda v:(means[str(v)],-v));return best,scores,means
def stability(feature_lists,rank_lists):
 out={'overall_full_selected':jaccard_sets(feature_lists),'overall_top20':jaccard_sets([x[:20] for x in rank_lists])}
 for b,p in PREFIXES.items():
  out[b]={'full_selected':jaccard_sets([[c for c in x if c.startswith(p)] for x in feature_lists]),'top20':jaccard_sets([[c for c in x if c.startswith(p)][:20] for x in rank_lists])}
 return out

def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--repo-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--seed',type=int,default=42);p.add_argument('--bootstrap-iterations',type=int,default=9999);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--chunk',type=int,default=5000);p.add_argument('--rsf-trees',type=int,default=500);p.add_argument('--max-steps',type=int,default=100);a=p.parse_args();start=time.perf_counter();root=a.input_dir.resolve();repo=a.repo_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 decisions={'exploratory':True,'cohort':'corrected primary-tumor-only','outer_folds':'identical to Table 7','variant_v':'top 200 per omics block by outer-training variance; clinical retained','variant_e':'top 2000 per omics block by outer-training variance, elastic-net Cox l1_ratio=0.5, nonzero coefficients refit by ridge CoxPH and RSF','elastic_alpha_grid':list(ALPHAS),'elastic_grid_deviation':'No repository elastic-net nested grid existed; grid fixed before results','ridge_grid':list(RIDGE_GRID),'inner_cv':'5-fold grouped/event-stratified, seed=4200+outer_fold','skipped_runtime':['Variant V DeepSurv']}
 base,y=load_cohort(root,'primary');folds=grouped_folds(y,5,a.seed);stored=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');stored=stored[(stored.cohort=='primary')&stored.configuration.isin(['Clinical-only','Full'])]
 task3=json.loads((out/'task3_fusion_strategies.json').read_text(encoding='utf8'));baseline_features=[f['baseline_fdr_selected_features'] for f in task3['folds']]
 fold_frames=[];records=[];scores={};vsets=[];vranks=[];esets=[];eranks=[]
 for fold,(tr,te) in enumerate(folds,1):
  xtr,xte,prep=prepare_fold(root,base,y,tr,te,a.meth_prefilter,a.knn_neighbors,a.chunk);yt,ye=y.iloc[tr],y.iloc[te];inner=inner_folds(yt,4200+fold,5);risks={};fr={'fold':fold,'preprocessing':prep,'configurations':{}};sd=stored[stored.fold==fold].set_index(['configuration','sample_key'])
  for cfg in ('Clinical-only','Full'):
   risks[cfg]={};fr['configurations'][cfg]={}
   for m in ('coxph','rsf'):
    r=sd.xs(cfg,level='configuration').loc[ye.index,m].to_numpy(float);risks[cfg][m]=r;fr['configurations'][cfg][m]={'cindex':ci(ye,r),'source':'stored Table 7 OOF'};scores.setdefault((cfg,m),[]).append(ci(ye,r))
  # Elastic-net comparison uses stored CoxPH references explicitly.
  risks['Clinical-only']['elastic_net_cox']=risks['Clinical-only']['coxph'];risks['Full']['elastic_net_cox']=risks['Full']['coxph']
  vsel,vrank=variance_select(xtr,200);vsets.append(vsel);vranks.append(sorted([c for b in PREFIXES for c in vrank[b][:200]],key=lambda c:float(xtr[c].var()),reverse=True));vt,ve=xtr[vsel],xte[vsel];cfg='Variant V — variance top200/block';risks[cfg]={};fr['configurations'][cfg]={'selected_features':vsel,'selected_counts':{b:sum(c.startswith(p) for c in vsel) for b,p in PREFIXES.items()}}
  bp,ins,means=tune_cox(vt,yt,inner,RIDGE_GRID,a.max_steps);cr,diag=cox_fit_predict(vt,yt,ve,bp,a.max_steps);rr=rsf_predict(vt,yt,ve,a.seed+fold,a.rsf_trees)
  for m,r in [('coxph',cr),('rsf',rr)]:risks[cfg][m]=r;fr['configurations'][cfg][m]={'cindex':ci(ye,r)};scores.setdefault((cfg,m),[]).append(ci(ye,r))
  fr['configurations'][cfg]['coxph'].update({'selected_penalizer':bp,'inner_scores':ins,'inner_means':means,'diagnostics':diag})
  esel0,erank0=variance_select(xtr,2000);et0=xtr[esel0];alpha,ains,ameans=tune_enet(et0,yt,inner);enet_r,coef,_=coxnet_fit(et0,yt,xte[esel0],alpha);nonzero=np.flatnonzero(np.abs(coef)>1e-12);esel=[esel0[i] for i in nonzero];order=[esel0[i] for i in np.argsort(np.abs(coef))[::-1] if abs(coef[i])>1e-12];esets.append(esel);eranks.append(order);cfg='Variant E — elastic-net nonzero';risks[cfg]={'elastic_net_cox':enet_r};fr['configurations'][cfg]={'variance_prefilter_counts':{b:min(2000,sum(c.startswith(p) for c in xtr)) for b,p in PREFIXES.items()},'selected_alpha':alpha,'inner_scores':ains,'inner_means':ameans,'selected_features':esel,'selected_nonzero':len(esel),'elastic_net_cox':{'cindex':ci(ye,enet_r)}};scores.setdefault((cfg,'elastic_net_cox'),[]).append(ci(ye,enet_r))
  if not esel:raise RuntimeError(f'Fold {fold}: elastic net selected zero features')
  bp,ins,means=tune_cox(xtr[esel],yt,inner,RIDGE_GRID,a.max_steps);cr,diag=cox_fit_predict(xtr[esel],yt,xte[esel],bp,a.max_steps);rr=rsf_predict(xtr[esel],yt,xte[esel],a.seed+fold,a.rsf_trees)
  for m,r in [('coxph',cr),('rsf',rr)]:risks[cfg][m]=r;fr['configurations'][cfg][m]={'cindex':ci(ye,r)};scores.setdefault((cfg,m),[]).append(ci(ye,r))
  fr['configurations'][cfg]['coxph'].update({'selected_penalizer':bp,'inner_scores':ins,'inner_means':means,'diagnostics':diag});records.append(fr);fold_frames.append({'fold':fold,'y':ye,'risks':risks});print('Task4 fold',fold,{c:{m:ci(ye,r) for m,r in z.items()} for c,z in risks.items()},flush=True)
 comparisons=[]
 for cfg in ('Variant V — variance top200/block','Variant E — elastic-net nonzero'):comparisons += [(f'{cfg} minus Clinical-only',cfg,'Clinical-only'),(f'{cfg} minus Full',cfg,'Full')]
 boot=paired_bootstrap(fold_frames,comparisons,a.bootstrap_iterations,a.seed);configs=list(dict.fromkeys(c for c,m in scores));payload={'task':'Phase 16 Task 4 alternative feature selection','status':'completed','decision_rules_fixed_before_results':decisions,'versions':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'sklearn':sklearn.__version__,'sksurv':sksurv.__version__,'lifelines':lifelines.__version__},'command':' '.join([sys.executable,*sys.argv]),'runtime_seconds':time.perf_counter()-start,'folds':records,'summary':{c:{m:summarize({m:v})[m] for (cc,m),v in scores.items() if cc==c} for c in configs},'bootstrap_comparisons':boot,'feature_stability':{'univariate_cox_baseline':stability(baseline_features,baseline_features),'variant_v':stability(vsets,vranks),'variant_e':stability(esets,eranks)}}
 (out/'task4_feature_selection_variants.json').write_text(json.dumps(payload,indent=2),encoding='utf8');lines=['# Phase 16 Task 4 — Alternative feature-selection methods','',f"Decision rules fixed before results: `{json.dumps(decisions)}`",'','## C-index','','| Configuration | Model | Mean | Sample SD | Fold values |','|---|---|---:|---:|---|']
 for c,models in payload['summary'].items():
  for m,z in models.items():lines.append(f"| {c} | {m} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_cindices'])}` |")
 lines+=['','## Bootstrap','','| Comparison | Model | Estimate | 95% CI | Raw p | BH q |','|---|---|---:|---|---:|---:|']
 for n,models in boot.items():
  for m,z in models.items():lines.append(f"| {n} | {m} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['raw_two_sided_p']} | {z['bh_fdr_q_within_family']} |")
 lines+=['','## Feature stability','',f"`{json.dumps(payload['feature_stability'])}`",'','## Reproducibility','',f"Command: `{payload['command']}`  ",f"Runtime seconds: `{payload['runtime_seconds']}`."];(out/'task4_feature_selection_variants.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
if __name__=='__main__':main()
