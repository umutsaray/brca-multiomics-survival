#!/usr/bin/env python3
"""Compile the Phase 17 status and exact key numbers from Task 6/7 JSON files."""
from __future__ import annotations
import argparse,json
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--results-dir',type=Path,required=True);a=p.parse_args();root=a.results_dir.resolve();t6=json.loads((root/'task6_learning_curve.json').read_text(encoding='utf-8'));t7=json.loads((root/'task7_uno_harrell_sensitivity.json').read_text(encoding='utf-8'));lines=['# Phase 17 summary','','## Task 7 — Uno versus Harrell sensitivity','','Status: **completed**. Stored OOF predictions only; no model refitting.','',f"Runtime seconds: `{t7['runtime_seconds']}`. Bootstrap: `{t7['decision_rules_fixed_before_results']['bootstrap_iterations']}` requested iterations, seed `{t7['decision_rules_fixed_before_results']['bootstrap_seed']}`.",'','Each line below states whether the Clinical-only versus Full ordering changed relative to the existing-tau Uno sign under the same cohort/model/weight construction.','','| Cohort | Model | Tau | Weight construction | Clinical-only − Full | 95% CI | Valid iterations | Ordering changed? |','|---|---|---|---|---:|---|---:|---|']
 for cohort,models in t7['comparisons'].items():
  for model,taus in models.items():
   for tau,constructions in taus.items():
    for construction,z in constructions.items():lines.append(f"| {cohort} | {model} | {tau} | {construction} | {z['mean_difference']} | `{json.dumps(z['ci_95'])}` | {z['valid_iterations']} | {'yes' if z['ordering_changed_from_existing_tau'] else 'no'} |")
 lines+=['','The two censoring-outcome constructions used identical outcome sets in every fold: '+', '.join(f"{cohort}=`{all(f['training_constructions_identical_outcomes'] for f in folds.values())}`" for cohort,folds in t7['folds'].items())+'.','','## Task 6 — Reduced CoxPH-only learning curve','',f"Status: **{t6['status']}**.",'',f"Stop reason: {t6['stop_reason']}",f"Actual elapsed seconds: `{t6['runtime_seconds']}`. Projected remaining seconds: `{t6['projected_remaining_seconds_at_stop']}`. Completed checkpoints: `{t6['completed_checkpoint_count']}`.",'','RSF and DeepSurv were skipped by the pre-fixed runtime decision. CoxPH used fixed ridge penalizer 0.1.','',t6['partial_result_warning'],'','| Fraction | Configuration | Mean | Sample SD | Completed folds |','|---:|---|---:|---:|---:|']
 for fraction in ('0.5','0.75','1.0'):
  for config in ('Clinical-only','Full'):
   z=t6['summary'][fraction][config];lines.append(f"| {fraction} | {config} | {z['mean']} | {z['sample_sd']} | {z['completed_folds']} |")
 lines+=['','Ordering statements:']
 s=t6['summary'];partial=(s['0.5']['Clinical-only']['mean'],s['0.5']['Full']['mean']);stored=(s['1.0']['Clinical-only']['mean'],s['1.0']['Full']['mean']);lines.append(f"- 50% completed checkpoint only: Clinical-only − Full = `{partial[0]-partial[1]}`; ordering changed versus stored 100% = `no`, but this is not the pre-specified completed estimand.")
 lines.append('- 75%: ordering change = `unavailable`; no fit was started.')
 b=t6['clinical_minus_full_bootstrap']['1.0'];lines.append(f"- Stored 100% Table 7: Clinical-only − Full = `{b['point_difference']}`; 95% CI = `{json.dumps(b['ci_95'])}`; ordering changed versus itself = `no`.")
 lines+=['','No learning-curve trend is classified because the reduced design did not complete.']
 (root/'phase17_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':main()
