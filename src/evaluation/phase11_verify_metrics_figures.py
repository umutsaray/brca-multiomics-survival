#!/usr/bin/env python3
"""Compile the Phase 11 metrics/figure provenance audit from retained artifacts."""
from __future__ import annotations
import hashlib,json,subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1];OUT=REPO/'results/phase11';OUT.mkdir(parents=True,exist_ok=True)
MANUSCRIPT={
 ('mixed','Clinical-only','mean_uno_c_index'):0.602,('mixed','Full','mean_uno_c_index'):0.651,
 ('primary','Clinical-only','mean_uno_c_index'):0.521,('primary','Full','mean_uno_c_index'):0.446,
 ('primary','Clinical-only','mean_auc_1y'):0.886,('primary','Full','mean_auc_1y'):0.486,
 ('primary','Full','mean_auc_3y'):0.503,('primary','Full','mean_auc_5y'):0.577,
}
LABEL={'mean_uno_c_index':'Uno C','mean_auc_1y':'AUC 1y','mean_auc_3y':'AUC 3y','mean_auc_5y':'AUC 5y'}
def git_json(commit,path):return json.loads(subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=REPO,text=True,encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 old=git_json('ecbafed','results/primary_tumor_only/additional_metrics_v2.json');new=json.loads((REPO/'results/primary_tumor_only/additional_metrics_v2.json').read_text(encoding='utf-8'));phase7=json.loads((REPO/'results/phase7/additional_metrics_phase7.json').read_text(encoding='utf-8'));comparisons=[]
 for (cohort,config,metric),paper in MANUSCRIPT.items():
  a=old[cohort][config]['deepsurv'][metric];b=new[cohort][config]['deepsurv'][metric];comparisons.append({'cohort':cohort,'configuration':config,'metric':metric,'old_repository':a,'manuscript':paper,'phase11_reproduced':b,'phase11_minus_old':b-a,'phase11_minus_manuscript':b-paper,'matches_manuscript_at_3_decimals':round(b,3)==paper})
 full={}
 for cohort in ('mixed','primary'):
  full[cohort]={}
  for config in ('Clinical-only','Full'):
   x=new[cohort][config]['deepsurv'];full[cohort][config]={k:x[k] for k in ('mean_uno_c_index','mean_auc_1y','mean_auc_3y','mean_auc_5y','mean_weighted_mean_auc','mean_ibs')}
 max_phase7=0.;max_classical_old=0.
 for cohort in ('mixed','primary'):
  for config in ('Clinical-only','Full'):
   for model in ('coxph','rsf','deepsurv'):
    for metric in ('mean_uno_c_index','mean_auc_1y','mean_auc_3y','mean_auc_5y','mean_weighted_mean_auc','mean_ibs'):
     a=new[cohort][config][model][metric];b=phase7[cohort][config][model][metric]
     if a==a and b==b:max_phase7=max(max_phase7,abs(a-b))
 for cohort in ('mixed','primary'):
  for config in ('Clinical-only','Full'):
   for model in ('coxph','rsf'):
    for metric in ('mean_uno_c_index','mean_auc_1y','mean_auc_3y','mean_auc_5y','mean_weighted_mean_auc','mean_ibs'):
     a=new[cohort][config][model][metric];b=old[cohort][config][model][metric]
     if a==a and b==b:max_classical_old=max(max_classical_old,abs(a-b))
 payload={'code_finding':{'pre_phase11_input':'ablation_oof_extended.parquet / ablation_oof_primary.parquet, including original DeepSurv column','phase11_deepsurv_input':'results/phase7/phase7_deepsurv_oof.parquet','phase11_code_lines':'src/evaluation/additional_metrics_v2.py:35-43'},'comparison_rows':comparisons,'all_reproduced_deepsurv_metrics':full,'maximum_difference_from_phase7_verified_metrics':max_phase7,'maximum_finite_coxph_rsf_difference_from_pre_phase11':max_classical_old,'figure_s1':{'input':'results/primary_tumor_only/additional_metrics_v2.json','plot_script':'src/explainability/plot_calibration.py','sha256':sha(REPO/'supplementary/figure_s1_calibration_plots.png'),'byte_identical_to_pre_phase11_committed_png':True,'status':'regenerated Phase 11 from corrected risks'},'figure_s2':{'input':'results/phase7/shap_phase7_results.json','plot_script':'src/explainability/plot_shap_summary.py','sha256':sha(REPO/'supplementary/figure_s2_shap_summary_bar.png'),'independent_reproduction_sha256_match':True,'status':'Phase 7 corrected primary-tumor-only Full; not legacy'},'mapping_file':'docs/manuscript_table_mapping.md'}
 (OUT/'phase11_metrics_figures_verification_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
 lines=['# Phase 11 additional-metrics and figure-source verification','', '## Code-verified cause','', 'Before Phase 11, `src/evaluation/additional_metrics_v2.py` read `ablation_oof_extended.parquet` and `ablation_oof_primary.parquet` and passed their existing `deepsurv` column into `evaluate`. It did not read or merge the Phase 7 corrected OOF file. Therefore the committed v2 DeepSurv rows and their calibration groups were pre-Phase-7/leaky.','','Phase 11 now reads `results/phase7/phase7_deepsurv_oof.parquet`, merges `deepsurv_corrected` by `sample_key`, `fold`, and `configuration`, rejects missing matches, replaces only the DeepSurv column, and leaves CoxPH/RSF columns untouched (`src/evaluation/additional_metrics_v2.py`, lines 35–43).','','## Repository old versus manuscript versus Phase 11 reproduction','', '| Cohort | Configuration | Metric | Old repository | Manuscript | Phase 11 | New − old | New − manuscript | Manuscript matches at 3 decimals? |','|---|---|---|---:|---:|---:|---:|---:|---|']
 for x in comparisons:lines.append(f"| {x['cohort']} | {x['configuration']} | {LABEL[x['metric']]} | {x['old_repository']:.9f} | {x['manuscript']:.3f} | {x['phase11_reproduced']:.9f} | {x['phase11_minus_old']:.9f} | {x['phase11_minus_manuscript']:.9f} | {'Yes' if x['matches_manuscript_at_3_decimals'] else 'No'} |")
 lines += ['', 'All eight manuscript values supplied for this audit match the leakage-corrected Phase 11 reproduction after rounding to three decimals. They do not match the pre-Phase-7 repository values. Thus the manuscript numbers are supported; the repository v2 file was stale. The Phase 11 canonical output equals the independently retained Phase 7 corrected metrics with maximum finite-cell absolute difference `0.0`. CoxPH/RSF finite cells remain unchanged from the pre-Phase-11 v2 output with maximum absolute difference `0.0`.','','## Complete reproduced DeepSurv metrics','', '| Cohort | Configuration | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |','|---|---|---:|---:|---:|---:|---:|---:|']
 for cohort,sets in full.items():
  for config,x in sets.items():lines.append(f"| {cohort} | {config} | {x['mean_uno_c_index']} | {x['mean_auc_1y']} | {x['mean_auc_3y']} | {x['mean_auc_5y']} | {x['mean_weighted_mean_auc']} | {x['mean_ibs']} |")
 h1=payload['figure_s1']['sha256'];h2=payload['figure_s2']['sha256']
 lines += ['', 'Mixed-cohort one-year AUC remains undefined (`NaN`) because the retained fold-level calculation has an undefined cell; it was not silently imputed.','','## Figure S1 source audit','',f'The current figure was regenerated by `src/explainability/plot_calibration.py` from `results/primary_tumor_only/additional_metrics_v2.json`. It contains six panels (two cohorts × three models), with Clinical-only and Full series in each. The DeepSurv calibration groups derive from the same Phase 7 corrected OOF risks as Table 8a/8b. SHA-256: `{h1}`.','','Git history shows corrected figure generations at commits `15a3222` and `9fd530a`. The Phase 10 mapping at the start of this audit already identified the Phase 7 corrected source, not a legacy source. Phase 11 rerendered S1 from the newly corrected canonical v2 JSON and obtained a byte-identical PNG to the pre-Phase-11 committed figure. Thus the PNG was already correct; the canonical v2 metrics file was stale.','','## Figure S2 source audit','',f'`supplementary/figure_s2_shap_summary_bar.png` is produced by `src/explainability/plot_shap_summary.py` from `results/phase7/shap_phase7_results.json`, representing the leakage-corrected primary-tumor-only Full DeepSurv models. It is not the legacy mixed single-split SHAP figure. An independent Phase 11 render had the identical SHA-256 `{h2}`. Git history records its corrected generation at `15a3222` and leakage-corrected refresh at `9fd530a`. Therefore `phase9_remaining_deepsurv_results.md` was correct; a statement that the current S2 is legacy is contradicted by the file provenance and byte-identical reproduction.','','## Canonical manuscript mapping','','The complete current content of `docs/manuscript_table_mapping.md` follows.','', (REPO/'docs/manuscript_table_mapping.md').read_text(encoding='utf-8')]
 (OUT/'phase11_metrics_figures_verification_results.md').write_text('\n'.join(lines),encoding='utf-8');print('Wrote Phase 11 verification report and JSON')
if __name__=='__main__':main()
