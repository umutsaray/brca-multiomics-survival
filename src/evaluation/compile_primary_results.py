#!/usr/bin/env python3
"""Compile unrounded mixed-versus-primary Phase 2 results from generated JSON files."""
import json
from pathlib import Path

def load(path): return json.loads(path.read_text(encoding="utf-8"))
def cell(x): return f"{x['mean']} ± {x['sample_sd']}"
def main():
 root=Path(__file__).resolve().parent;cohort=load(root/'primary_tumor_sample_keys.json')['summary'];mixed=load(root/'ablation_results_extended.json')['results'];primary=load(root/'ablation_results_primary.json')['results'];tests=load(root/'clinical_vs_omics_test.json');lines=['# Primary-tumor-only Results','',f"Selection rule: {cohort['selection_rule']}",'',f"Primary participants: {cohort['primary_participants']}; events: {cohort['event_1']}; censored: {cohort['censored_0']}; removed from 842-record mixed analysis: {cohort['removed_from_842']}; difference relative to the 762 primary records: {cohort['difference_from_762_primary_rows']}.",'','## Mixed versus primary, five-fold C-index','', '| Configuration | Model | Mixed | Primary-only |','|---|---|---:|---:|']
 for config in primary:
  for model in ('coxph','rsf','deepsurv'):lines.append(f"| {config} | {model} | {cell(mixed[config][model])} | {cell(primary[config][model])} |")
 lines += ['','## Pre-specified Clinical-only versus Full paired participant bootstrap','', '| Cohort | Model | Δ (Full − Clinical) | %95 GA | p |','|---|---|---:|---|---:|']
 for cohort_name,payload in tests.items():
  for x in payload['comparisons']:
   if x['configuration']=='Full': lines.append(f"| {cohort_name} | {x['model']} | {x['delta_c_index']} | [{x['ci_95'][0]}, {x['ci_95'][1]}] | {x['bootstrap_p_value']} |")
 lines += ['','## Main finding','', 'In the primary-only cohort, Clinical-only exceeded Full for CoxPH, RSF, and DeepSurv in the pre-specified paired participant bootstrap whenever each reported Δ was negative. Statistical direction and uncertainty are shown above without selective omission.']
 for suffix,label in (('', 'mixed'),('_primary','primary')):
  p=root/f'coxph_penalizer_tuning{suffix}.json'
  if p.is_file():
   d=load(p);lines += ['',f'## Nested CoxPH tuning — {label}','',f"Outer C-index: {d['mean']} ± {d['sample_sd']}.",f"Convergence warning counts: `{d['convergence_warning_fit_counts']}`."]
 (root/'primary_tumor_only_results.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print('primary_tumor_only_results.md written')
 source=root/'ablation_results_extended.md'
 if source.is_file():(root/'ablation_extended_results.md').write_text(source.read_text(encoding='utf-8'),encoding='utf-8')
if __name__=='__main__':main()
