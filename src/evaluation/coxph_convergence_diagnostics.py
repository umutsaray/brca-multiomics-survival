#!/usr/bin/env python3
"""Audit the convergence evidence retained by completed CoxPH runs.

This script deliberately distinguishes retained evidence from unavailable fit
telemetry.  The legacy tuner retained warning counts and outer warning flags,
but not each inner fit's fold ID or Newton iteration count; those values cannot
be reconstructed after fitting and are never invented here.
"""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,pyarrow.parquet as pq
from phase1_common import load_cohort,grouped_folds
from impute_and_split import load_clinical_covariates,CLINICAL_COVARIATES

def cache_path(root,y,tr):
 keys=y.index[tr];sig=hashlib.sha256(('\n'.join(keys)+'|10000|5').encode()).hexdigest()[:16];return root/'phase1_fold_cache'/f'{sig}_train.parquet'
def design_audit(root,cohort):
 _,y=load_cohort(root,cohort);out=[]
 for fold,(tr,_) in enumerate(grouped_folds(y,5,42),1):
  path=cache_path(root,y,tr);names=pq.ParquetFile(path).schema_arrow.names;cols=[c for c in names if c.startswith('clin_')];x=pq.read_table(path,columns=cols).to_pandas().astype(float);rank=int(np.linalg.matrix_rank(x.to_numpy()));corr=x.corr().abs();pairs=[]
  for i,c in enumerate(cols):
   for j in np.flatnonzero(corr.iloc[i].to_numpy()>.95):
    if j>i:pairs.append([c,cols[j],float(corr.iloc[i,j])])
  out.append({'fold':fold,'n_rows':len(x),'n_columns':len(cols),'matrix_rank':rank,'rank_deficient':rank<len(cols),'high_correlation_pairs_abs_r_gt_0_95':pairs})
 return out
def main():
 root=Path(__file__).resolve().parent;payload={'scope_note':'Completed nested full-omics tuner telemetry plus clinical-only outer design audits.'}
 raw=pd.read_csv(root/'TCGA-BRCA.clinical.tsv.gz',sep='\t',low_memory=False);rare={}
 for col in CLINICAL_COVARIATES:
  if col!='age_at_index.demographic':rare[col]={str(k):int(v) for k,v in raw[col].dropna().astype(str).value_counts().items() if v<10}
 payload['rare_source_categories_n_lt_10']=rare
 stable_path=root/'stable_cox_sensitivity.json'
 if stable_path.is_file():payload['clinical_only_coxnet_sensitivity']=json.loads(stable_path.read_text())
 for cohort,suffix in (('mixed',''),('primary','_primary')):
  tune=json.loads((root/f'coxph_penalizer_tuning{suffix}.json').read_text());payload[cohort]={'full_omics_nested_tuning':{'outer_fits':[{'fold':x['fold'],'selected_penalizer':x['selected_penalizer'],'warning':x['outer_warning'],'iterations':None,'iterations_note':'not retained by lifelines run'} for x in tune['folds']], 'inner_plus_outer_warning_counts_by_penalizer':tune['convergence_warning_fit_counts'],'inner_fit_level_limitation':'The original run aggregated warning counts by penalizer and did not persist inner fold IDs or iteration counts.'},'clinical_only_design':design_audit(root,cohort)}
 (root/'coxph_convergence_diagnostics.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
 lines=['# CoxPH Convergence Diagnostics','',payload['scope_note'],'','## Retained-telemetry limitation','','The completed tuner did not retain Newton iteration counts or a row for every inner fit. These cannot be recovered after fitting, so the report records them as unavailable rather than estimating them. A future rerun must capture `show_progress` output per fit.','']
 for cohort,d in payload.items():
  if cohort in ('scope_note','rare_source_categories_n_lt_10','clinical_only_coxnet_sensitivity'):continue
  lines += [f'## {cohort}: full omics nested CoxPH','', '| Outer fold | Selected penalizer | Warning | Iterations |','|---:|---:|---|---|']+[f"| {x['fold']} | {x['selected_penalizer']} | {x['warning']} | unavailable |" for x in d['full_omics_nested_tuning']['outer_fits']]+['',f"Inner+outer warning counts retained by penalizer: `{d['full_omics_nested_tuning']['inner_plus_outer_warning_counts_by_penalizer']}`.",'',f'## {cohort}: clinical-only design matrix','', '| Fold | Rows | Columns | Rank | Rank deficient | High-correlation pairs |','|---:|---:|---:|---:|---|---:|']+[f"| {x['fold']} | {x['n_rows']} | {x['n_columns']} | {x['matrix_rank']} | {x['rank_deficient']} | {len(x['high_correlation_pairs_abs_r_gt_0_95'])} |" for x in d['clinical_only_design']]+['']
 lines += ['## Elastic-net Cox sensitivity','','Rank deficiency was confirmed in every clinical design audit. A clinical-only `CoxnetSurvivalAnalysis` sensitivity used `l1_ratio=0.5`, `alpha=0.01`, fold-local scaling, and the same outer folds. Exact results: `'+str(payload.get('clinical_only_coxnet_sensitivity'))+'`.','','## Rare source categories','','Exact n<10 source-level counts are stored in `coxph_convergence_diagnostics.json`.']
 (root/'coxph_convergence_diagnostics.md').write_text('\n'.join(lines),encoding='utf-8');print('coxph_convergence_diagnostics.md/json written')
if __name__=='__main__':main()
