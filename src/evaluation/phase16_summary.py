#!/usr/bin/env python3
import json
from pathlib import Path
P=Path('results/phase16')
def load(n):return json.load(open(P/n,encoding='utf8'))
t0,t3,t4,t5=[load(x) for x in ('task0_her2_crosstab.json','task3_fusion_strategies.json','task4_feature_selection_variants.json','task5_shap_instability_decomposition.json')]
# Add the requested CoxPH design diagnostics to the Task 3 human-readable report.
t3_path=P/'task3_fusion_strategies.md';t3_text=t3_path.read_text(encoding='utf8')
if '## CoxPH design-matrix diagnostics' not in t3_text:
 D=['','## CoxPH design-matrix diagnostics','','| Fold | Configuration/block | Rows | Columns | Rank | Rank deficient | Convergence warnings |','|---:|---|---:|---:|---:|---|---|']
 for f in t3['folds']:
  for c,z in f['configurations'].items():
   if c=='Late fusion — equal Cox risks':
    for b,v in z['coxph']['blocks'].items():
     d=v['final_diagnostics'];D.append(f"| {f['fold']} | {c}: {b} | {d['n_rows']} | {d['n_columns']} | {d['matrix_rank']} | {d['rank_deficient']} | `{json.dumps(d['convergence_warnings'])}` |")
   elif isinstance(z,dict) and isinstance(z.get('coxph'),dict) and 'diagnostics' in z['coxph']:
    d=z['coxph']['diagnostics'];D.append(f"| {f['fold']} | {c} | {d['n_rows']} | {d['n_columns']} | {d['matrix_rank']} | {d['rank_deficient']} | `{json.dumps(d['convergence_warnings'])}` |")
 t3_path.write_text(t3_text+'\n'.join(D)+'\n',encoding='utf8')
# Keep raw bootstrap arrays in JSON while presenting a compact manuscript-facing report.
R=['# Phase 16 Task 4 — Alternative feature selection','',f"Status: {t4['status']}.",f"Stopped items: `{json.dumps(t4['stopped'])}`.",'','## C-index','','| Configuration | Model | Mean | Sample SD | Fold values |','|---|---|---:|---:|---|']
for c,ms in t4['summary'].items():
 for m,z in ms.items():R.append(f"| {c} | {m} | {z['mean']} | {z['sample_sd']} | `{json.dumps(z['fold_cindices'])}` |")
R+=['','## Fold-stratified paired bootstrap','','| Comparison | Model | Estimate | 95% CI | Raw p | BH q |','|---|---|---:|---|---:|---:|']
for n,ms in t4['bootstrap_comparisons'].items():
 for m,z in ms.items():R.append(f"| {n} | {m} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['raw_two_sided_p']} | {z['bh_fdr_q_within_family']} |")
R+=['','## Feature stability','',f"`{json.dumps(t4['feature_stability'])}`",'','## Reproducibility','',f"Command: `{t4['command']}`  ",f"Runtime seconds: `{t4['runtime_seconds']}`."]
(P/'task4_feature_selection_variants.md').write_text('\n'.join(R)+'\n',encoding='utf8')
L=['# Phase 16 summary','','Values are copied from Phase 16 JSON outputs without display rounding.','','## Task 0',f"- Status: {t0['status']}; HER2 conflicts={t0['conflict_participants']}; full tuple count sum={sum(x['n'] for x in t0['conflict_value_tuples'])}.",'','## Task 3',f"- Status: {t3['status']}; runtime_seconds={t3['runtime_seconds']}."]
for c,ms in t3['summary'].items():
 if c in ('Clinical-only','Full'):continue
 for m,z in ms.items():
  ref=t3['summary']['Clinical-only'][m]['mean'];delta=ref-z['mean'];L.append(f"- {c}, {m}: Clinical-only={ref}; new={z['mean']}; Clinical-only-new={delta}; ordering changed={'yes' if delta<0 else 'no'}.")
L+=['','## Task 4',f"- Status: {t4['status']}; runtime_seconds={t4['runtime_seconds']}; stopped={t4['stopped']}."]
for c,ms in t4['summary'].items():
 if c in ('Clinical-only','Full'):continue
 for m,z in ms.items():
  ref=t4['summary']['Clinical-only'][m]['mean'];delta=ref-z['mean'];L.append(f"- {c}, {m}: Clinical-only={ref}; new={z['mean']}; Clinical-only-new={delta}; ordering changed={'yes' if delta<0 else 'no'}.")
L+=['','## Task 5',f"- Status: {t5['status']}; runtime_seconds={t5['runtime_seconds']}; stored C-index reproduced exactly={all(x['cindex_difference']==0 for x in t5['folds'])}."]
for f in t5['folds']:L.append(f"- Fold {f['fold']}: within-fold top-20 Jaccard={f['within_fold_top20_jaccard']['mean_pairwise_jaccard']}; rank Spearman={f['within_fold_rank_spearman']['mean']}.")
L.append(f"- Selected in all five folds={t5['all_fold_overlap']['n_selected_in_all_five_folds']}; cross-fold common-feature rank Spearman={t5['all_fold_overlap']['cross_fold_rank_spearman']['mean']}.")
(P/'phase16_summary.md').write_text('\n'.join(L)+'\n',encoding='utf8');print(P/'phase16_summary.md')
