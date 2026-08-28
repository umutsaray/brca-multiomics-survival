# Manuscript table and result mapping

The manuscript itself was not present in the workspace. The labels below are those explicitly supplied in the Phase 3 request plus all numerical result families available locally. Any manuscript numbering not listed here requires author verification.

| Manuscript item | Producing script | Repository output | Status |
|---|---|---|---|
| Cohort construction/sample counts | `src/data_prep/build_primary_cohort.py` | `results/primary_tumor_only/primary_cohort_summary.md` | available |
| Table 2b — Clinical characteristics | `src/evaluation/clinical_characteristics.py` | `results/primary_tumor_only/clinical_characteristics.md` | available |
| Feature-count sensitivity / former Table 3b | `src/feature_selection/select_features_v2.py` | `results/mixed_cohort/select_features_v2_results.md`; `results/primary_tumor_only/select_features_v2_results_primary.md` | available |
| Table 4 — mixed-cohort ablation | `src/models/ablation_study.py --cohort mixed` | `results/mixed_cohort/ablation_results_extended.md` | available |
| Table 4b — primary-only ablation | `src/models/ablation_study.py --cohort primary` | `results/primary_tumor_only/ablation_results_primary.md`; `primary_tumor_only_results.md` | available |
| Table 5 — pre-specified Clinical-only vs Full test | `src/evaluation/clinical_vs_omics_test.py` | `results/primary_tumor_only/clinical_vs_omics_test.md` | available |
| Nested CoxPH penalizer result | `src/models/tune_coxph_penalizer.py` | `coxph_penalizer_tuning*.md/json` in both cohort result directories | available |
| Table 6b — convergence diagnostics | `src/evaluation/coxph_convergence_diagnostics.py` | `results/primary_tumor_only/coxph_convergence_diagnostics.md` | partial telemetry; iteration counts unavailable |
| OOF paired model comparison / former Table 7 | `src/evaluation/statistical_comparison_v2.py` | `statistical_comparison_v2*.md/json` | available |
| OOF score-scale diagnostic | `src/evaluation/oof_scale_diagnostics.py` | `results/primary_tumor_only/oof_scale_diagnostics.md` | training-mean centering unavailable |
| Table 9a/9b — Uno C, AUC, IBS, calibration | `src/evaluation/additional_metrics_v2.py` | `results/primary_tumor_only/additional_metrics_v2.md/json` | available; some mixed 1-year fold AUC values are NaN |
| Table 10 — primary-only SHAP | `src/explainability/explain_model_v2.py` | `results/primary_tumor_only/shap_v2_results.md/json` | available |
| Supplementary Table S1 | repository assembly from SHAP v2 JSON | `supplementary/table_s1_shap_full_ranking.csv` | available |
| Supplementary Table S2 | repository assembly from paired bootstrap JSON | `supplementary/table_s2_exploratory_clinical_vs_omics.csv` | available |
| Supplementary Table S3 | repository assembly from convergence JSON | `supplementary/table_s3_coxph_convergence_diagnostics.csv` | available |
| Supplementary Figure S1 | legacy `additional_metrics.py` | `supplementary/figure_s1_calibration_plots.png` | available; mixed legacy source, not primary v2 |
| Supplementary Figure S2 | legacy `explain_model.py` | `supplementary/figure_s2_shap_summary_bar.png` | available; superseded mixed single-split SHAP |
| Any other manuscript-numbered table | source manuscript not supplied | not mapped | author verification required |
