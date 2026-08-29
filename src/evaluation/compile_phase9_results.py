#!/usr/bin/env python3
"""Compile Phase 9 ablation summaries and refresh pooled exploratory comparisons."""
from __future__ import annotations
import argparse,json,sys
from datetime import date
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored
from statsmodels.stats.multitest import multipletests

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
sys.path.insert(0,str(HERE))

CONFIGS=['Clinical-only','RNA-seq only','Methylation only','miRNA only','Clinical + RNA-seq','Clinical + RNA-seq + Methylation','Clinical + RNA-seq + miRNA','Full','Clinical + Methylation','Clinical + miRNA']
REMAINING=[x for x in CONFIGS if x not in ('Clinical-only','Full')]

def ci(e,t,r):return float(concordance_index_censored(e.astype(bool),t.astype(float),r.astype(float))[0])

def paired_test_fast(data,iterations,seed):
 wide=data.pivot(index=['sample_key','fold','event','time'],columns='configuration',values=['coxph','rsf','deepsurv']);wide.columns=[f'{m}::{c}' for m,c in wide.columns];wide=wide.reset_index();comparisons=[];raw={}
 for model in ('coxph','rsf','deepsurv'):
  for config in CONFIGS:
   if config=='Clinical-only':continue
   a=f'{model}::{config}';b=f'{model}::Clinical-only';pair=wide.dropna(subset=[a,b]).reset_index(drop=True);participants=pair.sample_key.str[:12];groups=pair.groupby(participants,sort=True).indices;blocks=[np.asarray(x,dtype=int) for x in groups.values()];n=len(blocks);e=pair.event.to_numpy(bool);t=pair.time.to_numpy(float);ra=pair[a].to_numpy(float);rb=pair[b].to_numpy(float);point=ci(e,t,ra)-ci(e,t,rb);rng=np.random.default_rng(seed);arr=np.empty(iterations)
   for i in range(iterations):
    draw=rng.choice(n,size=n,replace=True);idx=np.concatenate([blocks[j] for j in draw]);arr[i]=ci(e[idx],t[idx],ra[idx])-ci(e[idx],t[idx],rb[idx])
   lo,hi=np.percentile(arr,[2.5,97.5]);pv=min(1.,2*min((np.count_nonzero(arr<=0)+1)/(iterations+1),(np.count_nonzero(arr>=0)+1)/(iterations+1)));role='pre-specified primary' if config=='Full' else 'exploratory';label=f'{model}: {config} - Clinical-only';raw[label]=arr.tolist();comparisons.append({'model':model,'configuration':config,'comparison_role':role,'delta_c_index':point,'ci_95':[float(lo),float(hi)],'contains_zero':bool(lo<=0<=hi),'bootstrap_p_value':float(pv),'n_pair_rows':len(pair),'n_pair_participants':n,'valid_bootstrap_iterations':iterations})
 exploratory=[i for i,r in enumerate(comparisons) if r['comparison_role']=='exploratory'];q=multipletests([comparisons[i]['bootstrap_p_value'] for i in exploratory],method='fdr_bh')[1]
 for i,v in zip(exploratory,q):comparisons[i]['bh_fdr_p_value']=float(v)
 for r in comparisons:
  if r['comparison_role']=='pre-specified primary':r['bh_fdr_p_value']=None
  r['interpretation']='anlamlı fark bulunamadı' if r['contains_zero'] else 'bootstrap %95 GA sıfırı içermiyor'
 return {'n_rows':len(wide),'n_participants':int(wide.sample_key.str[:12].nunique()),'iterations_requested':iterations,'seed':seed,'resampling_unit':'12-character TCGA participant ID','comparison_direction':'configuration minus Clinical-only','comparisons':comparisons,'raw_paired_differences':raw}

def corrected_grid(cohort,old,new7,new9):
 base=old[old.configuration.eq('Clinical-only')][['sample_key','fold','event','time']]
 grid=pd.concat([base.assign(configuration=c) for c in CONFIGS],ignore_index=True)
 classical=old[['sample_key','fold','event','time','configuration','coxph','rsf']]
 grid=grid.merge(classical,on=['sample_key','fold','event','time','configuration'],how='left',validate='one_to_one')
 dl=pd.concat([new7[(new7.cohort.eq(cohort)) & (new7.configuration.isin(['Clinical-only','Full']))],new9[new9.cohort.eq(cohort)]],ignore_index=True)
 grid=grid.merge(dl[['sample_key','fold','event','time','configuration','deepsurv']],on=['sample_key','fold','event','time','configuration'],how='left',validate='one_to_one')
 return grid

def fold_c(frame,col='deepsurv'):
 from sksurv.metrics import concordance_index_censored
 if frame.empty or frame[col].isna().all():return None
 z=frame.dropna(subset=[col]);return float(concordance_index_censored(z.event.astype(bool),z.time,z[col])[0])

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reuse-bootstrap',action='store_true',help='Reuse an existing Phase 9 bootstrap JSON while rebuilding reports.');args=parser.parse_args()
 out=REPO/'results/phase9';p9=json.loads((out/'phase9_deepsurv_results.json').read_text(encoding='utf-8'));new9=pd.read_parquet(out/'phase9_deepsurv_oof.parquet');new7=pd.read_parquet(REPO/'results/phase7/phase7_deepsurv_oof.parquet');old_test=json.loads((REPO/'results/primary_tumor_only/clinical_vs_omics_test.json').read_text(encoding='utf-8'))
 paths={'mixed':REPO/'results/mixed_cohort/ablation_oof_extended.parquet','primary':REPO/'results/primary_tumor_only/ablation_oof_primary.parquet'};grids={};summary={}
 lines=['# Phase 9 remaining DeepSurv ablation correction','',f"Requested fit attempts: {p9['fits_requested']}; completed: {p9['fits_completed']}; skipped because internal-train FDR selection yielded zero requested-block features: {p9['fits_skipped']}.",'','No fallback features and no old/leaky predictions were substituted for skipped folds. Thus every retained DeepSurv risk in the current ten-configuration grid comes from the leakage-corrected pipeline.','', '## Old versus corrected DeepSurv C-index','']
 for cohort,path in paths.items():
  old=pd.read_parquet(path);grid=corrected_grid(cohort,old,new7,new9);grids[cohort]=grid;summary[cohort]={};lines += [f'### {cohort}','','| Configuration | Fold | Old C | Corrected C | Corrected − old |','|---|---:|---:|---:|---:|']
  for config in REMAINING:
   olds=[];news=[];rows=[]
   for fold in range(1,6):
    oc=fold_c(old[(old.configuration.eq(config))&(old.fold.eq(fold))]);nc=fold_c(grid[(grid.configuration.eq(config))&(grid.fold.eq(fold))]);
    if oc is not None:olds.append(oc)
    if nc is not None:news.append(nc)
    rows.append({'fold':fold,'old_cindex':oc,'corrected_cindex':nc,'difference':None if oc is None or nc is None else nc-oc})
    lines.append(f"| {config} | {fold} | {'NA' if oc is None else f'{oc:.6f}'} | {'NA' if nc is None else f'{nc:.6f}'} | {'NA' if oc is None or nc is None else f'{nc-oc:.6f}'} |")
   summary[cohort][config]={'folds':rows,'old_n_folds':len(olds),'old_mean':float(np.mean(olds)) if olds else None,'old_sample_sd':float(np.std(olds,ddof=1)) if len(olds)>1 else None,'corrected_n_folds':len(news),'corrected_mean':float(np.mean(news)) if news else None,'corrected_sample_sd':float(np.std(news,ddof=1)) if len(news)>1 else None}
  lines += ['','Mean ± sample SD (means use only folds with an estimable model):','','| Configuration | Old | Corrected |','|---|---:|---:|']
  for c,s in summary[cohort].items():
   def cell(prefix):
    m=s[prefix+'_mean'];sd=s[prefix+'_sample_sd'];n=s[prefix+'_n_folds'];return 'NA' if m is None else f'{m:.6f} ± {sd:.6f} (n={n} folds)'
   lines.append(f'| {c} | {cell("old")} | {cell("corrected")} |')
  lines.append('')
 pd.concat([g.assign(cohort=c) for c,g in grids.items()],ignore_index=True).to_parquet(out/'ablation_oof_phase9.parquet',index=False)

 # Current Table 4/4b cells: unchanged CoxPH/RSF plus corrected DeepSurv for all ten configurations.
 ablation={};ab_md=['# Phase 9 current ablation tables','','CoxPH and RSF are unchanged retained estimates. Every reported DeepSurv value is leakage-corrected; Methylation-only means use only estimable folds and explicitly report n.','']
 for cohort,g in grids.items():
  ablation[cohort]={};ab_md += [f'## {cohort}','','| Configuration | CoxPH | RSF | DeepSurv corrected |','|---|---:|---:|---:|']
  for config in CONFIGS:
   ablation[cohort][config]={}
   cells=[]
   for model in ('coxph','rsf','deepsurv'):
    vals=[fold_c(g[(g.configuration.eq(config))&(g.fold.eq(f))],model) for f in range(1,6)];vals=[v for v in vals if v is not None];entry={'fold_cindices':vals,'n_folds':len(vals),'mean':float(np.mean(vals)) if vals else None,'sample_sd':float(np.std(vals,ddof=1)) if len(vals)>1 else None};ablation[cohort][config][model]=entry;cells.append('NA' if not vals else f"{entry['mean']:.6f} ± {entry['sample_sd']:.6f} (n={len(vals)})")
   ab_md.append(f"| {config} | {' | '.join(cells)} |")
  ab_md.append('')
 (out/'ablation_results_phase9.json').write_text(json.dumps(ablation,indent=2),encoding='utf-8');(out/'ablation_results_phase9.md').write_text('\n'.join(ab_md),encoding='utf-8')

 # Existing pooled-score method, intentionally retained for Supplementary Table S2.
 test_path=out/'clinical_vs_omics_phase9.json'
 if args.reuse_bootstrap and test_path.is_file():tests=json.loads(test_path.read_text(encoding='utf-8'))
 else:tests={c:paired_test_fast(g,1000,42) for c,g in grids.items()};test_path.write_text(json.dumps(tests,indent=2),encoding='utf-8')
 exploratory=[]
 for cohort,payload in tests.items():
  for r in payload['comparisons']:
   if r['comparison_role']=='exploratory':exploratory.append({'cohort':cohort,**r})
 pd.DataFrame(exploratory).drop(columns=['comparison_role']).to_csv(REPO/'supplementary/table_s2_exploratory_clinical_vs_omics.csv',index=False)
 test_lines=['# Phase 9 clinical-only versus omics pooled participant bootstrap','','This intentionally retains the prior pooled-score method for Supplementary Table S2. Inputs are Phase 9 leakage-corrected DeepSurv risks for all retained configurations and unchanged CoxPH/RSF risks. Iterations: 1,000; seed: 42.','']
 for cohort,payload in tests.items():
  test_lines += [f'## {cohort}','','| Role | Model | Configuration | Delta (configuration − Clinical) | 95% CI | p | BH-FDR p |','|---|---|---|---:|---:|---:|---:|']
  for r in payload['comparisons']:
   if r['comparison_role']!='exploratory':continue
   test_lines.append(f"| exploratory | {r['model']} | {r['configuration']} | {r['delta_c_index']:.6f} | [{r['ci_95'][0]:.6f}, {r['ci_95'][1]:.6f}] | {r['bootstrap_p_value']:.6f} | {r['bh_fdr_p_value']:.6f} |")
  test_lines.append('')
 (out/'clinical_vs_omics_phase9.md').write_text('\n'.join(test_lines),encoding='utf-8')

 # Direction/significance audit for the 16 changed exploratory DeepSurv inputs.
 changes=[]
 for cohort in ('mixed','primary'):
  oldmap={(r['model'],r['configuration']):r for r in old_test[cohort]['comparisons']};newmap={(r['model'],r['configuration']):r for r in tests[cohort]['comparisons']}
  for config in REMAINING:
   o=oldmap[('deepsurv',config)];n=newmap[('deepsurv',config)];changes.append({'cohort':cohort,'configuration':config,'old_delta':o['delta_c_index'],'new_delta':n['delta_c_index'],'old_ci':o['ci_95'],'new_ci':n['ci_95'],'old_p':o['bootstrap_p_value'],'new_p':n['bootstrap_p_value'],'old_bh_fdr_p':o['bh_fdr_p_value'],'new_bh_fdr_p':n['bh_fdr_p_value'],'direction_changed':bool(np.sign(o['delta_c_index'])!=np.sign(n['delta_c_index'])),'raw_ci_significance_changed':bool(o['contains_zero']!=n['contains_zero']),'bh_significance_changed':bool((o['bh_fdr_p_value']<.05)!=(n['bh_fdr_p_value']<.05))})
 lines += ['## Direction and significance audit','','Table 4/4b are descriptive and contain no inferential significance test; significance change is therefore not defined for those tables. Direction below refers to the pooled Supplementary Table S2 delta (configuration minus corrected Clinical-only). BH significance uses q < 0.05.','','| Cohort | Configuration | Old delta | New delta | Direction changed? | Raw-CI significance changed? | BH significance changed? |','|---|---|---:|---:|---|---|---|']
 for x in changes:lines.append(f"| {x['cohort']} | {x['configuration']} | {x['old_delta']:.6f} | {x['new_delta']:.6f} | {'Yes' if x['direction_changed'] else 'No'} | {'Yes' if x['raw_ci_significance_changed'] else 'No'} | {'Yes' if x['bh_significance_changed'] else 'No'} |")
 lines += ['','## SHAP / Table 10 impact check','','Table 10, Supplementary Table S1, and Figure S2 use only the primary-tumor-only Full DeepSurv models. Full was already retrained under the corrected pipeline in Phase 7 and is not one of the eight Phase 9 configurations. No SHAP recomputation is required; the Phase 7 SHAP outputs remain current.','','## Supplementary Table S2 impact','','The retained file contains 48 exploratory rows: 24 per cohort (8 configurations × 3 models). Exactly 16 rows have changed model inputs: the 8 DeepSurv exploratory comparisons in each of two cohorts. CoxPH and RSF inputs are unchanged, but their BH-adjusted q-values were recomputed jointly because the family includes the changed DeepSurv p-values. The complete updated table is `supplementary/table_s2_exploratory_clinical_vs_omics.csv`; full raw arrays are in `clinical_vs_omics_phase9.json`.','','## Method and seeds','','- Outer split: five-fold participant-grouped/event-stratified StratifiedGroupKFold, seed 42.','- Internal split: first held fold of seven-fold participant-grouped/event-stratified StratifiedGroupKFold, seed 4200 + outer fold.','- Feature selection, imputation, clinical encoding, and scaling: fitted on internal train only.','- DeepSurv training seed: 42 + outer fold; validation used only for early stopping; outer test predicted once.','- Network/optimizer/loss settings are identical to Phase 7.','- Supplementary Table S2: pooled participant bootstrap, 1,000 iterations, seed 42, BH correction over 24 exploratory comparisons within each cohort.','',f"Attempted {p9['fits_requested']} new fits; completed {p9['fits_completed']}; skipped {p9['fits_skipped']} zero-feature Methylation-only folds. Together with the 20 Phase 7 fits, 100 corrected fit attempts cover all ten configurations; 95 models were estimable and fitted.",f"Generation date: {date.today().isoformat()}. Git provenance: the commit containing this report; its immutable hash is reported in the delivery message because a commit cannot embed its own hash."]
 payload={'generation_date':date.today().isoformat(),'git_provenance':'commit containing this report; hash reported in delivery','method':p9['method'],'fits_requested':p9['fits_requested'],'fits_completed':p9['fits_completed'],'fits_skipped':p9['fits_skipped'],'ablation_comparison':summary,'supplementary_table_s2':{'exploratory_rows':48,'rows_with_changed_deepsurv_inputs':16,'changes':changes},'shap_impact':{'requires_recomputation':False,'reason':'SHAP uses primary Full, already corrected in Phase 7'}}
 (out/'phase9_remaining_deepsurv_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');(out/'phase9_remaining_deepsurv_results.md').write_text('\n'.join(lines),encoding='utf-8');print('Wrote Phase 9 reports, corrected grid, and updated Supplementary Table S2')

if __name__=='__main__':main()
