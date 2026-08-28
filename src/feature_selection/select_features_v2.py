#!/usr/bin/env python3
"""BH-FDR reporting and leakage-safe feature-count sensitivity analysis."""
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
from phase1_common import load_cohort,grouped_folds,prepare_fold,rank_fold,select_from_rankings,fit_models,add_model_args

def args():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--cohort',choices=['mixed','primary'],default='mixed');p.add_argument('--counts',type=int,nargs='+',default=[50,100,200,500]);p.add_argument('--folds',type=int,default=5);p.add_argument('--random-state',type=int,default=42);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--chunk-size',type=int,default=512);p.add_argument('--methylation-chunksize',type=int,default=5000);add_model_args(p);return p.parse_args()
def main():
 a=args();root=a.input_dir.resolve();base,y=load_cohort(root,a.cohort);suffix='_primary' if a.cohort=='primary' else '';results={str(n):{'coxph':[],'rsf':[],'deepsurv':[]} for n in a.counts};fdr=[];rna_sets=[];oof=[];started=time.perf_counter()
 for fold,(tr,te) in enumerate(grouped_folds(y,a.folds,a.random_state),1):
  fs=time.perf_counter();xtr,xte,prep=prepare_fold(root,base,y,tr,te,a.meth_prefilter,a.knn_neighbors,a.methylation_chunksize);ytr,yte=y.iloc[tr],y.iloc[te];ranks=rank_fold(xtr,ytr,a.chunk_size);fd={}
  for b,r in ranks.items():
   q=np.full(len(r),np.nan);ok=r.p_value.notna();q[ok]=multipletests(r.loc[ok,'p_value'],method='fdr_bh')[1];fd[b]=int(np.sum(q<.05));print(f'Fold {fold} {b}: FDR<0.05 = {fd[b]} / {int(ok.sum())}')
  fdr.append(fd)
  for n in a.counts:
   selected,stats=select_from_rankings(ranks,xtr,n,fdr=False);models=fit_models(xtr[selected],ytr,xte[selected],yte,a,a.random_state+fold*100+n)
   for m in results[str(n)]: results[str(n)][m].append(models[m]['cindex'])
   print(f'Fold={fold} N={n} '+ ' '.join(f'{m}={models[m]["cindex"]!r}' for m in models))
   if n==200:
    rna_sets.append(set(c for c in selected if c.startswith('rna_')))
    for j,idx in enumerate(te): oof.append({'sample_key':y.index[idx],'event':int(y.iloc[idx].event),'time':float(y.iloc[idx].time),**{m:float(models[m]['risk'][j]) for m in models},'fold':fold})
  print(f'Fold {fold} elapsed_seconds={time.perf_counter()-fs!r}')
 summary={}
 for n,v in results.items(): summary[n]={m:{'fold_cindices':s,'mean':float(np.mean(s)),'sample_sd':float(np.std(s,ddof=1))} for m,s in v.items()}
 pairs=[len(rna_sets[i]&rna_sets[j])/len(rna_sets[i]|rna_sets[j]) for i in range(len(rna_sets)) for j in range(i+1,len(rna_sets))]
 stability={'five_fold_intersection':len(set.intersection(*rna_sets)),'pairwise_jaccard':pairs,'mean_pairwise_jaccard':float(np.mean(pairs))}
 payload={'n_samples':len(y),'participant_count':int(pd.Index(y.index).str[:12].nunique()),'fdr_counts_by_fold':fdr,'sensitivity':summary,'rna_stability_top200':stability,'elapsed_seconds':time.perf_counter()-started}
 payload['cohort']=a.cohort;(root/f'select_features_v2_results{suffix}.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');pd.DataFrame(oof).to_parquet(root/f'phase1_oof_predictions{suffix}.parquet',index=False)
 lines=['# FDR and Feature-count Sensitivity','','Outer folds are event-stratified and participant-grouped. All preprocessing is fold-local.','','| N/block | CoxPH | RSF | DeepSurv |','|---:|---:|---:|---:|']
 for n,s in summary.items(): lines.append('| '+n+' | '+' | '.join(f"{s[m]['mean']} ± {s[m]['sample_sd']}" for m in ['coxph','rsf','deepsurv'])+' |')
 lines+=['',f"RNA five-fold intersection (top-200): {stability['five_fold_intersection']}",f"RNA mean pairwise Jaccard: {stability['mean_pairwise_jaccard']}",'',f"Elapsed seconds: {payload['elapsed_seconds']}"]
 (root/f'select_features_v2_results{suffix}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
if __name__=='__main__':main()
