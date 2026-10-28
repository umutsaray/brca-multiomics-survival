#!/usr/bin/env python3
"""Phase 16 Task 3: leakage-safe alternative fusion strategies."""
from __future__ import annotations
import argparse,json,platform,sys,time,warnings
from pathlib import Path
import lifelines,numpy as np,pandas as pd,sklearn,sksurv
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
HERE=Path(__file__).resolve().parent
for sub in ('common','data_prep','preprocessing','feature_selection','models','evaluation'):sys.path.insert(0,str(HERE.parents[0]/sub))
from phase1_common import load_cohort,grouped_folds,prepare_fold,rank_fold,select_from_rankings
from phase16_common import RIDGE_GRID,ci,cox_fit_predict,inner_folds,tune_cox,rsf_predict,summarize,paired_bootstrap

BLOCKS={'clinical':'clin_','rna':'rna_','mirna':'mirna_','methylation':'meth_'}
K_VALUES=(5,10,20)

def block_oof(x,y,cols,folds,grid,max_steps):
 p,scores,means=tune_cox(x[cols],y,folds,grid,max_steps);oof=np.empty(len(y))
 for tr,va in folds:oof[va]=cox_fit_predict(x.iloc[tr][cols],y.iloc[tr],x.iloc[va][cols],p,max_steps)[0]
 return oof,p,scores,means
def standardize_train_test(a,b):
 mu=float(np.mean(a));sd=float(np.std(a));sd=sd if sd>0 else 1.;return (a-mu)/sd,(b-mu)/sd,{'mean':mu,'sd':sd}
def pca_frame(xtr,xte,selected,k):
 parts_tr=[];parts_te=[];meta={}
 clinical=[c for c in xtr if c.startswith('clin_')];parts_tr.append(xtr[clinical].reset_index(drop=True));parts_te.append(xte[clinical].reset_index(drop=True))
 for block,prefix in [('rna','rna_'),('mirna','mirna_'),('methylation','meth_')]:
  cols=[c for c in selected if c.startswith(prefix)];sc=StandardScaler();a=sc.fit_transform(xtr[cols]);b=sc.transform(xte[cols]);n=min(k,len(cols),len(xtr));pc=PCA(n_components=n,random_state=42);u=pc.fit_transform(a);v=pc.transform(b);names=[f'{block}_pc{i+1}' for i in range(n)];parts_tr.append(pd.DataFrame(u,columns=names));parts_te.append(pd.DataFrame(v,columns=names));meta[block]={'selected_features':len(cols),'components':n,'explained_variance_ratio':pc.explained_variance_ratio_.tolist()}
 a=pd.concat(parts_tr,axis=1);b=pd.concat(parts_te,axis=1);a.index=xtr.index;b.index=xte.index;return a,b,meta

def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--repo-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--seed',type=int,default=42);p.add_argument('--bootstrap-iterations',type=int,default=9999);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--chunk',type=int,default=5000);p.add_argument('--top-n',type=int,default=200);p.add_argument('--rsf-trees',type=int,default=500);p.add_argument('--max-steps',type=int,default=100);a=p.parse_args();start=time.perf_counter();root=a.input_dir.resolve();repo=a.repo_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 decisions={'exploratory':True,'cohort':'corrected primary-tumor-only','outer_folds':'identical to Table 7','outer_seed':a.seed,'selection':'outer-training-only univariate Cox, BH-FDR < 0.05, cap 200/block','ridge_grid':list(RIDGE_GRID),'inner_cv':'5-fold participant-grouped/event-stratified; seed=4200+outer_fold','late_equal':'mean of standardized clinical/RNA/miRNA/methylation Cox risk scores','late_stack':'ridge CoxPH meta-model on four standardized inner-OOF risks','pca_k_values':list(K_VALUES),'pca_fit_scope':'outer-training only, separately per omics block','skipped_runtime':['RSF-based block late fusion','PCA DeepSurv']}
 base,y=load_cohort(root,'primary');folds=grouped_folds(y,5,a.seed);stored=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');stored=stored[(stored.cohort=='primary')&stored.configuration.isin(['Clinical-only','Full'])]
 fold_records=[];fold_frames=[];score_lists={};selected_by_fold=[]
 for fold,(tr,te) in enumerate(folds,1):
  xtr,xte,prep=prepare_fold(root,base,y,tr,te,a.meth_prefilter,a.knn_neighbors,a.chunk);yt,ye=y.iloc[tr],y.iloc[te];rank=rank_fold(xtr,yt);baseline_selected,selmeta=select_from_rankings(rank,xtr,a.top_n,True);baseline_selected=[c for c in baseline_selected if c in xtr];selected=[]
  for block in ('rna','mirna','methylation'):
   q=rank[block];selected+=q[q.p_value.notna()].head(min(a.top_n,int(q.p_value.notna().sum()))).feature.tolist()
  selected += [c for c in xtr if c.startswith('clin_')];selected_by_fold.append(selected);inner=inner_folds(yt,4200+fold,5);risks={};fr={'fold':fold,'preprocessing':prep,'selection':selmeta,'selected_features':selected,'baseline_fdr_selected_features':baseline_selected,'configurations':{}}
  sd=stored[stored.fold==fold].set_index(['configuration','sample_key'])
  for cfg in ('Clinical-only','Full'):
   risks[cfg]={};fr['configurations'][cfg]={}
   for model in ('coxph','rsf'):
    r=sd.xs(cfg,level='configuration').loc[ye.index,model].to_numpy(float);risks[cfg][model]=r;fr['configurations'][cfg][model]={'cindex':ci(ye,r),'source':'stored Table 7 OOF'};score_lists.setdefault((cfg,model),[]).append(ci(ye,r))
  block_train={};block_test={};block_meta={}
  for block,prefix in BLOCKS.items():
   cols=[c for c in selected if c.startswith(prefix)];oof,best,ins,means=block_oof(xtr,yt,cols,inner,RIDGE_GRID,a.max_steps);test,diag=cox_fit_predict(xtr[cols],yt,xte[cols],best,a.max_steps);u,v,scale=standardize_train_test(oof,test);block_train[block]=u;block_test[block]=v;block_meta[block]={'n_features':len(cols),'selected_penalizer':best,'inner_scores':ins,'inner_means':means,'risk_standardization':scale,'final_diagnostics':diag}
  A=pd.DataFrame(block_train,index=xtr.index);B=pd.DataFrame(block_test,index=xte.index);eq=B.mean(axis=1).to_numpy();cfg='Late fusion — equal Cox risks';risks[cfg]={'coxph':eq};fr['configurations'][cfg]={'coxph':{'cindex':ci(ye,eq),'blocks':block_meta}};score_lists.setdefault((cfg,'coxph'),[]).append(ci(ye,eq))
  mp,mins,mmeans=tune_cox(A,yt,inner,RIDGE_GRID,a.max_steps);stack,mdiag=cox_fit_predict(A,yt,B,mp,a.max_steps);cfg='Late fusion — stacked Cox risks';risks[cfg]={'coxph':stack};fr['configurations'][cfg]={'coxph':{'cindex':ci(ye,stack),'selected_penalizer':mp,'inner_scores':mins,'inner_means':mmeans,'diagnostics':mdiag}};score_lists.setdefault((cfg,'coxph'),[]).append(ci(ye,stack))
  for k in K_VALUES:
   pt,pe,pmeta=pca_frame(xtr,xte,selected,k);cfg=f'PCA early fusion k={k}';risks[cfg]={};fr['configurations'][cfg]={}
   bp,ins,means=tune_cox(pt,yt,inner,RIDGE_GRID,a.max_steps);cr,diag=cox_fit_predict(pt,yt,pe,bp,a.max_steps);rr=rsf_predict(pt,yt,pe,a.seed+fold,a.rsf_trees)
   for model,r in [('coxph',cr),('rsf',rr)]:risks[cfg][model]=r;fr['configurations'][cfg][model]={'cindex':ci(ye,r)};score_lists.setdefault((cfg,model),[]).append(ci(ye,r))
   fr['configurations'][cfg]['coxph'].update({'selected_penalizer':bp,'inner_scores':ins,'inner_means':means,'diagnostics':diag});fr['configurations'][cfg]['pca']=pmeta
  fold_records.append(fr);fold_frames.append({'fold':fold,'y':ye,'risks':risks});print('Task3 fold',fold,{c:{m:ci(ye,r) for m,r in z.items()} for c,z in risks.items()},flush=True)
 configs=[c for c,_ in score_lists];configs=list(dict.fromkeys(configs));comparisons=[]
 for cfg in configs:
  if cfg in ('Clinical-only','Full'):continue
  comparisons.extend([(f'Clinical-only minus {cfg}','Clinical-only',cfg),(f'{cfg} minus Full',cfg,'Full')])
 boot=paired_bootstrap(fold_frames,comparisons,a.bootstrap_iterations,a.seed);payload={'task':'Phase 16 Task 3 alternative fusion strategies','status':'completed','decision_rules_fixed_before_results':decisions,'versions':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'sklearn':sklearn.__version__,'sksurv':sksurv.__version__,'lifelines':lifelines.__version__},'command':' '.join([sys.executable,*sys.argv]),'runtime_seconds':time.perf_counter()-start,'folds':fold_records,'summary':{c:{m:summarize({m:v})[m] for (cc,m),v in score_lists.items() if cc==c} for c in configs},'bootstrap_comparisons':boot}
 (out/'task3_fusion_strategies.json').write_text(json.dumps(payload,indent=2),encoding='utf8')
 lines=['# Phase 16 Task 3 — Alternative fusion strategies','',f"Decision rules fixed before results: `{json.dumps(decisions)}`",'','## C-index','','| Configuration | Model | Mean | Sample SD | Fold values |','|---|---|---:|---:|---|']
 for c,models in payload['summary'].items():
  for m,z in models.items():lines.append(f"| {c} | {m} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_cindices'])}` |")
 lines+=['','## Fold-stratified paired bootstrap','','| Comparison | Model | Estimate | 95% CI | Raw p | BH q |','|---|---|---:|---|---:|---:|']
 for n,models in boot.items():
  for m,z in models.items():lines.append(f"| {n} | {m} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['raw_two_sided_p']} | {z['bh_fdr_q_within_family']} |")
 lines+=['','## Reproducibility','',f"Command: `{payload['command']}`  ",f"Runtime seconds: `{payload['runtime_seconds']}`."];(out/'task3_fusion_strategies.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
if __name__=='__main__':main()
