# Manuscript table and figure mapping

This is the Phase 11 canonical mapping. “Current” identifies the source that should populate the manuscript; superseded pooled, leaky, or historical single-split outputs remain available only for audit.

| Manuscript item | Producing script | Canonical repository output | Status / provenance |
|---|---|---|---|
| **Table 1.** Comparative literature overview | Author literature compilation; present-study rows: `src/evaluation/compile_primary_results.py` | `results/primary_tumor_only/primary_tumor_only_results.md` | Current; literature-derived except present-study rows |
| **Table 2.** Raw omics layers and cohort intersection | `src/data_prep/prepare_dataset.py` | `results/mixed_cohort/data_summary.md` | Current |
| **Table 2b.** Clinical characteristics | `src/evaluation/clinical_characteristics.py` | `results/primary_tumor_only/clinical_characteristics.md`; `.json` | Current; both cohorts |
| **Table 3.** Original single-split univariate Cox selection summary | `src/feature_selection/select_features.py` | `results/mixed_cohort/selected_features.json` | Historical selection summary retained for traceability |
| **Table 3b.** Fold-local feature-count sensitivity | `src/feature_selection/select_features_v2.py` | `results/mixed_cohort/select_features_v2_results.*`; `results/primary_tumor_only/select_features_v2_results_primary.*` | Current sensitivity analysis |
| **Table 4.** Mixed-cohort ten-configuration ablation | Phase 7: `src/models/phase7_deepsurv.py`; Phase 9: `src/models/phase9_deepsurv.py`; compilation: `src/evaluation/compile_phase9_results.py` | `results/phase9/ablation_results_phase9.md`; `.json` | Phase 9 current for all ten configurations; zero-feature Methylation-only folds are explicit NA |
| **Table 4b.** Primary-tumor-only ten-configuration ablation | Same as Table 4 | `results/phase9/ablation_results_phase9.md`; `.json` | Phase 9 current; leakage-corrected DeepSurv throughout |
| **Table 5.** Pre-specified Clinical-only versus Full comparison | `src/evaluation/foldstratified_test.py`; `src/evaluation/compile_phase7_results.py` | `results/phase7/foldstratified_test_results.md`; `.json` | Current: participant-level fold-stratified bootstrap, 9,999 iterations |
| **Table 6.** CoxPH convergence and design-matrix diagnostics | `src/evaluation/coxph_convergence_diagnostics.py` | `results/primary_tumor_only/coxph_convergence_diagnostics.md`; `.json` | Current; formerly numbered Table 6b |
| **Table 7.** Leakage-corrected Full DeepSurv versus CoxPH/RSF | `src/evaluation/table8_fix.py` | `results/phase8/table8_fix_results.md`; `.json` | Current Phase 8 fold-stratified result; formerly Table 8 |
| **Table 8a.** Mixed-cohort Uno C, time-dependent AUC, IBS | `src/evaluation/additional_metrics_v2.py` | `results/primary_tumor_only/additional_metrics_v2.md`; `.json` | Phase 11 canonical; DeepSurv uses Phase 7 corrected OOF; CoxPH/RSF unchanged |
| **Table 8b.** Primary-tumor-only Uno C, time-dependent AUC, IBS | `src/evaluation/additional_metrics_v2.py` | `results/primary_tumor_only/additional_metrics_v2.md`; `.json` | Phase 11 canonical; DeepSurv uses Phase 7 corrected OOF; CoxPH/RSF unchanged |
| **Table 9.** Primary-tumor-only Full DeepSurv SHAP ranking | `src/models/phase7_deepsurv.py`; `src/evaluation/compile_phase7_results.py` | `results/phase7/shap_phase7_results.md`; `.json` | Current leakage-corrected result; formerly Table 10 |
| **Supplementary Table S1.** Full 20-feature SHAP ranking | `src/models/phase7_deepsurv.py`; `src/evaluation/compile_phase7_results.py` | `supplementary/table_s1_shap_full_ranking.csv`; `results/phase7/shap_phase7_results.json` | Current Phase 7 ranking |
| **Supplementary Table S2.** Exploratory Clinical-only versus omics comparisons | `src/evaluation/phase10_exploratory_foldstratified.py` | `supplementary/table_s2_exploratory_clinical_vs_omics.csv`; `results/phase10/phase10_exploratory_foldstratified_results.*` | Current Phase 10 fold-stratified method; supersedes Phase 9 pooled-score version |
| **Supplementary Table S3.** Fold-level CoxPH diagnostics | `src/evaluation/coxph_convergence_diagnostics.py` | `supplementary/table_s3_coxph_convergence_diagnostics.csv` | Current diagnostic supplement |
| **Supplementary Table S4.** Historical uncorrected single-split model comparison | `src/models/train_baselines.py`; `src/models/train_explainable_dl.py`; `src/evaluation/cross_validate.py` | `results/mixed_cohort/baseline_results.json`; `model_comparison.md`; `cv_results.*` | Historical/uncorrected; formerly main-text Table 6 |
| **Supplementary Table S5.** Historical Wilcoxon/single-split bootstrap comparison | `src/evaluation/statistical_comparison.py`; `src/evaluation/bootstrap_ci.py` | `results/mixed_cohort/statistical_comparison.*`; `bootstrap_ci_results.*` | Historical/uncorrected; formerly main-text Table 7 |
| **Supplementary Figure S1.** Five-year calibration, both cohorts/configurations/models | Metrics: `src/evaluation/additional_metrics_v2.py`; plot: `src/explainability/plot_calibration.py` | `supplementary/figure_s1_calibration_plots.png`; provenance readme; input: `results/primary_tumor_only/additional_metrics_v2.json` | Phase 11 regenerated from corrected DeepSurv OOF; six panels, two configurations per panel |
| **Supplementary Figure S2.** Primary-tumor-only Full DeepSurv SHAP bar plot | SHAP: `src/models/phase7_deepsurv.py`; plot: `src/explainability/plot_shap_summary.py` | `supplementary/figure_s2_shap_summary_bar.png`; provenance readme; input: `results/phase7/shap_phase7_results.json` | Current Phase 7 leakage-corrected model; SHA-256 reproduction verified in Phase 11; not legacy single-split SHAP |

## Superseded sources

- `results/primary_tumor_only/additional_metrics_v2.*` before Phase 11 used DeepSurv risks from the original ablation OOF files and is superseded by the same canonical paths regenerated in Phase 11.
- `results/phase9/clinical_vs_omics_phase9.*` uses pooled scores and is superseded by Phase 10 for Supplementary Table S2.
- Original mixed/primary ablation JSON files remain the CoxPH/RSF source and historical DeepSurv audit source; current DeepSurv ablation cells are in Phase 9.
