#!/usr/bin/env python3
"""Leakage-safe five-fold omics-block ablation study."""
import argparse,json,time
from pathlib import Path
import numpy as np
from phase1_common import load_cohort,grouped_folds,prepare_fold,rank_fold,select_from_rankings,fit_models,add_model_args
COMBOS={'Clinical-only':['clin_'],'RNA-seq only':['rna_'],'Methylation only':['meth_'],'miRNA only':['mirna_'],'Clinical + RNA-seq':['clin_','rna_'],'Clinical + RNA-seq + Methylation':['clin_','rna_','meth_'],'Clinical + RNA-seq + miRNA':['clin_','rna_','mirna_'],'Full':['clin_','rna_','meth_','mirna_'],'Clinical + Methylation':['clin_','meth_'],'Clinical + miRNA':['clin_','mirna_']}
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--cohort',choices=['mixed','primary'],default='mixed');p.add_argument('--top-n',type=int,default=200);p.add_argument('--random-state',type=int,default=42);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--methylation-chunksize',type=int,default=5000);add_model_args(p);a=p.parse_args();root=a.input_dir.resolve();base,y=load_cohort(root,a.cohort);suffix='_primary' if a.cohort=='primary' else '_extended';res={k:{m:[] for m in ['coxph','rsf','deepsurv']} for k in COMBOS};oof=[];feature_counts=[];start=time.perf_counter()
 for fold,(tr,te) in enumerate(grouped_folds(y,5,a.random_state),1):
  xtr,xte,_=prepare_fold(root,base,y,tr,te,a.meth_prefilter,a.knn_neighbors,a.methylation_chunksize);ytr,yte=y.iloc[tr],y.iloc[te];rank=rank_fold(xtr,ytr);sel,stats=select_from_rankings(rank,xtr,a.top_n,True);print(f'Fold {fold} FDR selection: {stats}')
  for name,prefixes in COMBOS.items():
   for block,prefix in {'rna':'rna_','mirna':'mirna_','methylation':'meth_','clinical':'clin_'}.items():
    if prefix in prefixes:
     count=len([c for c in sel if c.startswith(prefix)])
     feature_counts.append({'cohort':a.cohort,'config':name,'fold':fold,'omics_block':block,'n_features_selected':count})
   cols=[c for c in sel if c.startswith(tuple(prefixes))];
   if not cols: print(f'Fold {fold} {name}: no FDR-qualified features; skipped');continue
   models=fit_models(xtr[cols],ytr,xte[cols],yte,a,a.random_state+fold);print(f'Fold {fold} {name}: '+str({m:z['cindex'] for m,z in models.items()}))
   for j,key in enumerate(yte.index): oof.append({'sample_key':key,'fold':fold,'configuration':name,'event':int(yte.iloc[j].event),'time':float(yte.iloc[j].time),**{m:float(models[m]['risk'][j]) for m in models}})
   for m in res[name]:res[name][m].append(models[m]['cindex'])
 summary={k:{m:{'fold_cindices':v,'mean':float(np.mean(v)) if v else None,'sample_sd':float(np.std(v,ddof=1)) if len(v)>1 else None} for m,v in z.items()} for k,z in res.items()};payload={'cohort':a.cohort,'results':summary,'elapsed_seconds':time.perf_counter()-start};(root/f'ablation_results{suffix}.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');lines=['# Omics Ablation Study','','| Combination | CoxPH | RSF | DeepSurv |','|---|---:|---:|---:|'];
 for k,z in summary.items():lines.append('| '+k+' | '+' | '.join('NA' if z[m]['mean'] is None else f"{z[m]['mean']} ± {z[m]['sample_sd']}" for m in ['coxph','rsf','deepsurv'])+' |')
 (root/f'ablation_results{suffix}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');pd=__import__('pandas');pd.DataFrame(oof).to_parquet(root/f'ablation_oof{suffix}.parquet',index=False);fc=pd.DataFrame(feature_counts);fc.to_csv(root/f'ablation_feature_counts{suffix}.csv',index=False)
 audit=['# Ablation feature-count audit','',f'- Cohort: {a.cohort}',f'- Outer folds: 5-fold StratifiedGroupKFold, participant-grouped, event-stratified, shuffle seed {a.random_state}.',f'- Selection: fold-local univariate Cox, BH-FDR < 0.05, capped at {a.top_n} per omics block.','- No FDR-threshold relaxation or feature fallback exists. A zero-feature single-block configuration is skipped; in a combined configuration, the empty block contributes no columns and the remaining blocks are fitted.','', '| Configuration | Fold | Block | Selected features |','|---|---:|---|---:|']
 audit += [f"| {r.config} | {r.fold} | {r.omics_block} | {r.n_features_selected} |" for r in fc.itertuples(index=False)]
 meth=fc[fc.omics_block.eq('methylation')];audit += ['',f"Methylation-containing fold/configuration rows with zero selected probes: {int(meth.n_features_selected.eq(0).sum())}.",f"Methylation-containing fold/configuration rows with fewer than five selected probes: {int(meth.n_features_selected.lt(5).sum())}.",'']
 (root/f'ablation_feature_count_audit{suffix}.md').write_text('\n'.join(audit),encoding='utf-8')
if __name__=='__main__':main()
