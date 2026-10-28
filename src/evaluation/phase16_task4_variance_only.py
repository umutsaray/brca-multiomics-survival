#!/usr/bin/env python3
"""Complete Phase 16 Task 4 Variant V; record Variant E runtime stop."""
from __future__ import annotations
import argparse,json,platform,sys,time
from pathlib import Path
import lifelines,numpy as np,pandas as pd,sklearn,sksurv
HERE=Path(__file__).resolve().parent
for sub in ('common','data_prep','preprocessing','feature_selection','models','evaluation'):sys.path.insert(0,str(HERE.parents[0]/sub))
from phase1_common import load_cohort,grouped_folds,prepare_fold
from phase16_common import RIDGE_GRID,ci,cox_fit_predict,inner_folds,tune_cox,rsf_predict,summarize,paired_bootstrap,jaccard_sets
P={'rna':'rna_','mirna':'mirna_','methylation':'meth_'}
def stability(sets,ranks):
 d={'overall_full_selected':jaccard_sets(sets),'overall_top20':jaccard_sets([x[:20] for x in ranks])}
 for b,p in P.items():d[b]={'full_selected':jaccard_sets([[c for c in x if c.startswith(p)] for x in sets]),'top20':jaccard_sets([[c for c in x if c.startswith(p)][:20] for x in ranks])}
 return d
def main():
 q=argparse.ArgumentParser();q.add_argument('--input-dir',type=Path,required=True);q.add_argument('--repo-dir',type=Path,required=True);q.add_argument('--output-dir',type=Path,required=True);q.add_argument('--seed',type=int,default=42);q.add_argument('--bootstrap-iterations',type=int,default=9999);q.add_argument('--rsf-trees',type=int,default=500);a=q.parse_args();st=time.perf_counter();root=a.input_dir.resolve();repo=a.repo_dir.resolve();out=a.output_dir.resolve();task3=json.loads((out/'task3_fusion_strategies.json').read_text());base,y=load_cohort(root,'primary');stored=pd.read_parquet(repo/'results/phase9/ablation_oof_phase9.parquet');stored=stored[(stored.cohort=='primary')&stored.configuration.isin(['Clinical-only','Full'])];frames=[];folds=[];scores={};sets=[];ranks=[]
 for fold,(tr,te) in enumerate(grouped_folds(y,5,a.seed),1):
  xtr,xte,prep=prepare_fold(root,base,y,tr,te,10000,5,5000);yt,ye=y.iloc[tr],y.iloc[te];sel=[];rank=[]
  for b,p in P.items():s=xtr[[c for c in xtr if c.startswith(p)]].var().sort_values(ascending=False);sel+=s.index[:200].tolist();rank+=s.index.tolist()
  sel += [c for c in xtr if c.startswith('clin_')];sets.append(sel);ranks.append(rank);risks={};sd=stored[stored.fold==fold].set_index(['configuration','sample_key']);fr={'fold':fold,'preprocessing':prep,'selected_features':sel,'selected_counts':{b:sum(c.startswith(p) for c in sel) for b,p in P.items()},'configurations':{}}
  for cfg in ('Clinical-only','Full'):
   risks[cfg]={}
   for m in ('coxph','rsf'):
    r=sd.xs(cfg,level='configuration').loc[ye.index,m].to_numpy();risks[cfg][m]=r;scores.setdefault((cfg,m),[]).append(ci(ye,r))
  cfg='Variant V — variance top200/block';inner=inner_folds(yt,4200+fold);bp,ins,means=tune_cox(xtr[sel],yt,inner,RIDGE_GRID,100);cr,diag=cox_fit_predict(xtr[sel],yt,xte[sel],bp,100);rr=rsf_predict(xtr[sel],yt,xte[sel],a.seed+fold,a.rsf_trees);risks[cfg]={'coxph':cr,'rsf':rr};fr['configurations'][cfg]={'coxph':{'cindex':ci(ye,cr),'selected_penalizer':bp,'inner_scores':ins,'inner_means':means,'diagnostics':diag},'rsf':{'cindex':ci(ye,rr)}}
  for m,r in risks[cfg].items():scores.setdefault((cfg,m),[]).append(ci(ye,r))
  folds.append(fr);frames.append({'fold':fold,'y':ye,'risks':risks});print('Task4V fold',fold,fr['configurations'][cfg],flush=True)
 comps=[('Variant V minus Clinical-only','Variant V — variance top200/block','Clinical-only'),('Variant V minus Full','Variant V — variance top200/block','Full')];boot=paired_bootstrap(frames,comps,a.bootstrap_iterations,a.seed);configs=list(dict.fromkeys(c for c,m in scores));baseline=[f['baseline_fdr_selected_features'] for f in task3['folds']];payload={'task':'Phase 16 Task 4 alternative feature selection','status':'partially completed','stopped':{'Variant E':'Stopped before any outer-fold result was produced; projected runtime exceeded the approximately 4-hour notification/feasibility threshold after the first 6000-feature elastic-net inner path did not finish in approximately 12 minutes. No substitute analysis was used.','Variant V DeepSurv':'Skipped for runtime as pre-specified.'},'decision_rules_fixed_before_results':{'variant_v':'top 200/block by outer-training variance','ridge_grid':list(RIDGE_GRID),'outer_folds':'Table 7','seed':a.seed},'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__,'sksurv':sksurv.__version__,'lifelines':lifelines.__version__},'command':' '.join([sys.executable,*sys.argv]),'runtime_seconds':time.perf_counter()-st,'folds':folds,'summary':{c:{m:summarize({m:v})[m] for (cc,m),v in scores.items() if cc==c} for c in configs},'bootstrap_comparisons':boot,'feature_stability':{'univariate_cox_baseline':stability(baseline,baseline),'variant_v':stability(sets,ranks)}};(out/'task4_feature_selection_variants.json').write_text(json.dumps(payload,indent=2));lines=['# Phase 16 Task 4 — Alternative feature selection','','Status: partially completed. Variant E stopped for projected runtime >4 h; no substitute was used. Variant V DeepSurv skipped for runtime.','','## C-index','','| Configuration | Model | Mean | Sample SD | Fold values |','|---|---|---:|---:|---|']
 for c,ms in payload['summary'].items():
  for m,z in ms.items():lines.append(f"| {c} | {m} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_cindices'])}` |")
 lines+=['','## Bootstrap','',f"`{json.dumps(boot)}`",'','## Feature stability','',f"`{json.dumps(payload['feature_stability'])}`",'','## Reproducibility','',f"Command: `{payload['command']}`  ",f"Runtime seconds: `{payload['runtime_seconds']}`."];(out/'task4_feature_selection_variants.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
