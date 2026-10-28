#!/usr/bin/env python3
"""Phase 16 Task 5: decompose feature-selection and SHAP attribution instability."""
from __future__ import annotations
import argparse,json,platform,sys,time
from pathlib import Path
import numpy as np,pandas as pd,scipy,shap,torch
from scipy.stats import rankdata,spearmanr
HERE=Path(__file__).resolve().parent
for sub in ('common','data_prep','preprocessing','feature_selection','models','evaluation','explainability'):sys.path.insert(0,str(HERE.parents[0]/sub))
from phase1_common import load_cohort,grouped_folds
from phase7_deepsurv import fit_network,inner_split,shap_fold
from phase16_common import jaccard_sets
PREFIXES={'rna':'rna_','mirna':'mirna_','methylation':'meth_','clinical':'clin_'}
def mean_pairwise_spearman(rank_maps,features):
 vals=[]
 for i in range(len(rank_maps)):
  for j in range(i+1,len(rank_maps)):
   r=float(spearmanr([rank_maps[i][x] for x in features],[rank_maps[j][x] for x in features]).statistic) if len(features)>1 else None;vals.append({'pair':[i+1,j+1],'spearman':r})
 good=[x['spearman'] for x in vals if x['spearman'] is not None and np.isfinite(x['spearman'])];return {'pairs':vals,'mean':float(np.mean(good)) if good else None}
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--repo-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--seed',type=int,default=42);p.add_argument('--shap-seeds',type=int,default=10);p.add_argument('--background',type=int,default=30);p.add_argument('--explain',type=int,default=30);p.add_argument('--epochs',type=int,default=300);p.add_argument('--patience',type=int,default=30);p.add_argument('--lr',type=float,default=.001);p.add_argument('--weight-decay',type=float,default=.0001);p.add_argument('--dropout',type=float,default=.3);a=p.parse_args();torch.set_num_threads(1);start=time.perf_counter();root=a.input_dir.resolve();repo=a.repo_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 decisions={'exploratory':True,'model':'existing leakage-corrected primary Full DeepSurv','fixed_feature_set_within_fold':True,'shap_repetitions_per_fold':a.shap_seeds,'shap_seed_schema':'10000 + 100*outer_fold + repetition','background_n':a.background,'explained_n':a.explain,'model_seed_schema':'42 + outer_fold','internal_seed_schema':'4200 + outer_fold'}
 old=json.loads((repo/'results/phase7/shap_phase7_results.json').read_text(encoding='utf8'));base,y=load_cohort(root,'primary');folds=grouped_folds(y,5,a.seed);fold_out=[];selection_sets=[];fold_mean_importance=[]
 for fold,(otr,ote) in enumerate(folds,1):
  fit_idx,val_idx=inner_split(y,otr,4200+fold);yf,yv,yt=y.iloc[fit_idx],y.iloc[val_idx],y.iloc[ote];stem=f'primary_fold{fold}_innerseed{4200+fold}';cache=root/'phase7_inner_cache';xf=pd.read_parquet(cache/f'{stem}_fit.parquet');xv=pd.read_parquet(cache/f'{stem}_validation.parquet');xt=pd.read_parquet(cache/f'{stem}_test.parquet');selected_set=set(old['folds'][fold-1]['importance']);selected=[c for c in xf.columns if c in selected_set];xf,xv,xt=xf[selected],xv[selected],xt[selected];selection_sets.append(selected)
  fit=fit_network(xf,yf,xv,yv,xt,yt,a.seed+fold,a.epochs,a.patience,a.lr,a.weight_decay,a.dropout);runs=[]
  for rep in range(a.shap_seeds):
   seed=10000+100*fold+rep;s=shap_fold(fit,seed,a.background,a.explain);s['seed']=seed;runs.append(s)
  rank_maps=[]
  for s in runs:
   ordered=sorted(s['importance'],key=s['importance'].get,reverse=True);rank_maps.append({name:i+1 for i,name in enumerate(ordered)})
  jac=jaccard_sets([s['top20'] for s in runs]);sp=mean_pairwise_spearman(rank_maps,selected);meanimp={f:float(np.mean([s['importance'][f] for s in runs])) for f in selected};fold_mean_importance.append(meanimp)
  stored=float(old['folds'][fold-1].get('cindex',json.loads((repo/'results/phase7/phase7_deepsurv_results.json').read_text(encoding='utf8'))['cohorts']['primary']['folds'][fold-1]['models']['Full']['cindex']))
  fold_out.append({'fold':fold,'selected_features':selected,'n_selected':len(selected),'stored_cindex':stored,'reproduced_cindex':fit['cindex'],'cindex_difference':fit['cindex']-stored,'best_epoch':fit['best_epoch'],'stored_best_epoch':old['folds'][fold-1]['best_epoch'],'runs':runs,'within_fold_top20_jaccard':jac,'within_fold_rank_spearman':sp});print('Task5 fold',fold,fit['cindex'],jac['mean_pairwise_jaccard'],sp['mean'],flush=True)
 selection={'overall':jaccard_sets(selection_sets)}
 for b,pref in PREFIXES.items():selection[b]=jaccard_sets([[x for x in s if x.startswith(pref)] for s in selection_sets])
 common=sorted(set.intersection(*map(set,selection_sets)));fold_rank=[]
 for imp in fold_mean_importance:
  ordered=sorted(common,key=lambda f:imp[f],reverse=True);fold_rank.append({f:i+1 for i,f in enumerate(ordered)})
 overlap={'n_selected_in_all_five_folds':len(common),'features':common,'fewer_than_20':len(common)<20,'cross_fold_rank_spearman':mean_pairwise_spearman(fold_rank,common)}
 payload={'task':'Phase 16 Task 5 SHAP instability decomposition','status':'completed','decision_rules_fixed_before_results':decisions,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'shap':shap.__version__,'torch':torch.__version__},'command':' '.join([sys.executable,*sys.argv]),'runtime_seconds':time.perf_counter()-start,'folds':fold_out,'selection_induced_instability':selection,'all_fold_overlap':overlap,'stored_original_cross_fold_shap_top20_jaccard':old['mean_pairwise_jaccard']}
 (out/'task5_shap_instability_decomposition.json').write_text(json.dumps(payload,indent=2),encoding='utf8');lines=['# Phase 16 Task 5 — SHAP instability decomposition','',f"Decision rules fixed before results: `{json.dumps(decisions)}`",'','## Within-fold attribution noise','','| Fold | Stored C-index | Reproduced C-index | Difference | Mean top-20 Jaccard | Mean rank Spearman |','|---:|---:|---:|---:|---:|---:|']
 for z in fold_out:lines.append(f"| {z['fold']} | {z['stored_cindex']} | {z['reproduced_cindex']} | {z['cindex_difference']} | {z['within_fold_top20_jaccard']['mean_pairwise_jaccard']} | {z['within_fold_rank_spearman']['mean']} |")
 lines+=['','## Selection-induced instability','',f"`{json.dumps(selection)}`",'','## Features selected in all folds','',f"Count: `{len(common)}`; fewer than 20: `{len(common)<20}`.  ",f"Cross-fold rank consistency: `{json.dumps(overlap['cross_fold_rank_spearman'])}`.",'','## Reproducibility','',f"Command: `{payload['command']}`  ",f"Runtime seconds: `{payload['runtime_seconds']}`."];(out/'task5_shap_instability_decomposition.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
if __name__=='__main__':main()
