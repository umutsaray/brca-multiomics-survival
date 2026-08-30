# Manuscript table and figure mapping

This Phase 12 mapping follows sequential MDPI IJMS main-text numbering (`Table 1` through `Table 13`). The change is limited to manuscript labels and cross-references; producing code, retained outputs, analyses, and numerical results are unchanged.

## Main-text tables

| Previous number | Current number | Content / manuscript role | Producing script | Canonical repository output | Status / provenance |
|---|---|---|---|---|---|
| Table 1 | **Table 1** | Comparative summary of multi-omics BRCA studies in the literature | Author literature compilation; present-study rows: `src/evaluation/compile_primary_results.py` | `results/primary_tumor_only/primary_tumor_only_results.md` | Current; literature-derived except present-study rows |
| Table 2 | **Table 2** | Raw omics-layer characteristics and cohort intersection | `src/data_prep/prepare_dataset.py` | `results/mixed_cohort/data_summary.md` | Current |
| Table 2b | **Table 3** | Clinical characteristics of the mixed and primary-tumor-only cohorts | `src/evaluation/clinical_characteristics.py` | `results/primary_tumor_only/clinical_characteristics.md`; `.json` | Current; both cohorts |
| Table 3 | **Table 4** | Univariate Cox feature-selection summary for the original single split | `src/feature_selection/select_features.py` | `results/mixed_cohort/selected_features.json` | Historical selection summary retained for traceability |
| Table 3b | **Table 5** | Sensitivity of participant-grouped five-fold CV C-index to the number of retained features | `src/feature_selection/select_features_v2.py` | `results/mixed_cohort/select_features_v2_results.*`; `results/primary_tumor_only/select_features_v2_results_primary.*` | Current sensitivity analysis |
| Table 4 | **Table 6** | Ten-configuration ablation study, mixed cohort | Phase 7: `src/models/phase7_deepsurv.py`; Phase 9: `src/models/phase9_deepsurv.py`; compilation: `src/evaluation/compile_phase9_results.py` | `results/phase9/ablation_results_phase9.md`; `.json` | Phase 9 current for all configurations; zero-feature Methylation-only folds are explicit NA |
| Table 4b | **Table 7** | Ten-configuration ablation study, corrected primary-tumor-only cohort | Same as Table 6 | `results/phase9/ablation_results_phase9.md`; `.json` | Phase 9 current; leakage-corrected DeepSurv throughout |
| Table 5 | **Table 8** | Pre-specified primary comparison: fold-stratified paired-bootstrap C-index difference between Clinical-only and Full | `src/evaluation/foldstratified_test.py`; `src/evaluation/compile_phase7_results.py` | `results/phase7/foldstratified_test_results.md`; `.json` | Current; participant-level fold-stratified bootstrap, 9,999 iterations |
| Table 6 (formerly Table 6b) | **Table 9** | Clinical-only CoxPH design-matrix and convergence diagnostics | `src/evaluation/coxph_convergence_diagnostics.py` | `results/primary_tumor_only/coxph_convergence_diagnostics.md`; `.json` | Current |
| Table 7 (formerly Table 8) | **Table 10** | Leakage-corrected Full DeepSurv versus CoxPH and RSF, fold-stratified paired bootstrap | `src/evaluation/table8_fix.py` | `results/phase8/table8_fix_results.md`; `.json` | Current Phase 8 result |
| Table 8a (formerly Table 9a) | **Table 11** | Additional metrics—Uno’s C, time-dependent AUC, and IBS—mixed cohort | `src/evaluation/additional_metrics_v2.py` | `results/primary_tumor_only/additional_metrics_v2.md`; `.json` | Phase 11 canonical; corrected DeepSurv OOF, unchanged CoxPH/RSF |
| Table 8b (formerly Table 9b) | **Table 12** | Additional metrics—Uno’s C, time-dependent AUC, and IBS—primary-tumor-only cohort | `src/evaluation/additional_metrics_v2.py` | `results/primary_tumor_only/additional_metrics_v2.md`; `.json` | Phase 11 canonical; corrected DeepSurv OOF, unchanged CoxPH/RSF |
| Table 9 (formerly Table 10) | **Table 13** | Selected top-ranked SHAP features from the primary-tumor-only Full DeepSurv model | `src/models/phase7_deepsurv.py`; `src/evaluation/compile_phase7_results.py` | `results/phase7/shap_phase7_results.md`; `.json` | Current leakage-corrected result |

## Supplementary tables

Supplementary numbering is unchanged.

| Current number | Content | Producing script | Canonical repository output | Status / relationship to main text |
|---|---|---|---|---|
| **Table S1** | Full 20-feature SHAP ranking | `src/models/phase7_deepsurv.py`; `src/evaluation/compile_phase7_results.py` | `supplementary/table_s1_shap_full_ranking.csv`; `results/phase7/shap_phase7_results.json` | Current Phase 7 ranking; supports Table 13 |
| **Table S2** | Exploratory Clinical-only versus omics comparisons | `src/evaluation/phase10_exploratory_foldstratified.py` | `supplementary/table_s2_exploratory_clinical_vs_omics.csv`; `results/phase10/phase10_exploratory_foldstratified_results.*` | Current Phase 10 fold-stratified method |
| **Table S3** | Fold-level CoxPH diagnostics | `src/evaluation/coxph_convergence_diagnostics.py` | `supplementary/table_s3_coxph_convergence_diagnostics.csv` | Current diagnostic supplement; supports Table 9 |
| **Table S4** | Historical uncorrected single-split model comparison | Historical sources: `src/models/train_baselines.py`, `src/models/train_explainable_dl.py`, `src/evaluation/bootstrap_ci.py`; Phase 13 compilation: `src/evaluation/compile_phase13_supplementary.py` | `supplementary/table_s4_historical_single_split_evaluation.csv` | Historical/uncorrected; superseded by Table 6 and Table 7; compiled for submission in Phase 13 |
| **Table S5** | Historical Wilcoxon/single-split bootstrap comparison | Historical sources: `src/evaluation/statistical_comparison.py`, `src/evaluation/bootstrap_ci.py`; Phase 13 compilation: `src/evaluation/compile_phase13_supplementary.py` | `supplementary/table_s5_historical_wilcoxon_bootstrap_comparison.csv` | Historical/uncorrected; superseded by Table 10; compiled for submission in Phase 13 |

## Supplementary figures

Figure numbering and figure content are unchanged.

| Current number | Content | Producing script and input | Canonical repository output | Status / table cross-reference |
|---|---|---|---|---|
| **Figure S1** | Five-year calibration for both cohorts, both configurations, and all three models | Metrics: `src/evaluation/additional_metrics_v2.py`; plot: `src/explainability/plot_calibration.py`; input: `results/primary_tumor_only/additional_metrics_v2.json` | `supplementary/figure_s1_calibration_plots.png`; provenance readme | Phase 11 verified; accompanies the additional metrics in Table 11 and Table 12 |
| **Figure S2** | Primary-tumor-only Full DeepSurv SHAP summary bar plot | SHAP: `src/models/phase7_deepsurv.py`; plot: `src/explainability/plot_shap_summary.py`; input: `results/phase7/shap_phase7_results.json` | `supplementary/figure_s2_shap_summary_bar.png`; provenance readme | Phase 11 verified; accompanies the selected SHAP features in Table 13 |

## Superseded analysis sources retained for audit

- Pre-Phase-11 `results/primary_tumor_only/additional_metrics_v2.*` used the original DeepSurv risks; the canonical files at the same paths were regenerated in Phase 11 and now supply Table 11 and Table 12.
- `results/phase9/clinical_vs_omics_phase9.*` uses pooled scores and is superseded by Phase 10 for Table S2.
- Original mixed/primary ablation JSON files remain the unchanged CoxPH/RSF source and historical DeepSurv audit source; current DeepSurv ablation cells for Table 6 and Table 7 are in Phase 9.
