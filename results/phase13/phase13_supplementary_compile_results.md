# Phase 13 supplementary compilation results

## Scope

This phase only compiled retained historical values and packaged existing supplementary artefacts. No model was fitted and no statistical result was recomputed.

## Sources and exact extracted fields

- `results/mixed_cohort/baseline_results.json`
  - `models.coxph.test_c_index` = `0.7064765300059418`
  - `models.random_survival_forest.test_c_index` = `0.714795008912656`
- `results/mixed_cohort/model_comparison.md`
  - Retained single-split comparison row confirms DeepSurv test C-index; the unrounded value below is taken from the paired-bootstrap JSON.
- `results/mixed_cohort/cv_results.json`
  - Read for requested historical provenance verification (`folds` = `5`); no CV value is inserted into Table S4 or S5.
- `results/mixed_cohort/bootstrap_ci_results.json`
  - `models.CoxPH.test_c_index` = `0.7064765300059418`
  - `models.CoxPH.ci_95_percentile_lower/upper` = `[0.5146235268728914, 0.8622735382074371]`
  - `models.Random Survival Forest.test_c_index` = `0.714795008912656`
  - `models.Random Survival Forest.ci_95_percentile_lower/upper` = `[0.5598458138821741, 0.8279464235842526]`
  - `models.DeepSurv.test_c_index` = `0.7456922162804516`
  - `models.DeepSurv.ci_95_percentile_lower/upper` = `[0.6187689806562355, 0.8447662489579691]`
- `results/mixed_cohort/statistical_comparison.json`
  - `comparisons.DeepSurv vs CoxPH.wilcoxon_two_sided_p_value` = `0.3125`
  - `comparisons.DeepSurv vs CoxPH.bootstrap_mean_difference` = `0.03934835766237502`
  - `comparisons.DeepSurv vs CoxPH.bootstrap_95_ci_lower/upper` = `[-0.05316103010958706, 0.1663211835726834]`
  - `comparisons.DeepSurv vs CoxPH.significant_by_percentile_ci` = `false`
  - `comparisons.DeepSurv vs Random Survival Forest.wilcoxon_two_sided_p_value` = `0.125`
  - `comparisons.DeepSurv vs Random Survival Forest.bootstrap_mean_difference` = `0.03131173145862762`
  - `comparisons.DeepSurv vs Random Survival Forest.bootstrap_95_ci_lower/upper` = `[-0.06476056379129216, 0.1259146529730082]`
  - `comparisons.DeepSurv vs Random Survival Forest.significant_by_percentile_ci` = `false`

The retained bootstrap metadata are `requested_iterations = 1000` and `random_state = 42`. The retained statistical-comparison metadata are `bootstrap_iterations = 1000`, `paired_resampling = true`, and `cv_fold_count = 5`.

## Table S4 (full CSV content)

| Model | C-index | 95% CI | Validation scheme |
|---|---|---|---|
| CoxPH | 0.7064765300059418 | [0.5146235268728914, 0.8622735382074371] | Original, non-participant-grouped single train/test split |
| Random Survival Forest | 0.714795008912656 | [0.5598458138821741, 0.8279464235842526] | Original, non-participant-grouped single train/test split |
| DeepSurv | 0.7456922162804516 | [0.6187689806562355, 0.8447662489579691] | Original, non-participant-grouped single train/test split |

## Table S5 (full CSV content)

| Comparison | Wilcoxon p-value | Bootstrap mean difference | 95% CI | Significant |
|---|---|---|---|---|
| DeepSurv vs CoxPH | 0.3125 | 0.03934835766237502 | [-0.05316103010958706, 0.1663211835726834] | No |
| DeepSurv vs Random Survival Forest | 0.125 | 0.03131173145862762 | [-0.06476056379129216, 0.1259146529730082] | No |

## Submission archive

Archive: `supplementary/supplementary_materials.zip`

| File | Uncompressed bytes | Compressed bytes |
|---|---:|---:|
| `table_s1_shap_full_ranking.csv` | 988 | 458 |
| `annotated_shap_features.md` | 4649 | 1369 |
| `annotated_shap_features.json` | 7372 | 1313 |
| `table_s2_exploratory_clinical_vs_omics.csv` | 7727 | 2727 |
| `table_s3_coxph_convergence_diagnostics.csv` | 936 | 245 |
| `table_s4_historical_single_split_evaluation.csv` | 442 | 247 |
| `table_s5_historical_wilcoxon_bootstrap_comparison.csv` | 275 | 207 |
| `figure_s1_calibration_plots.png` | 512708 | 454472 |
| `figure_s2_shap_summary_bar.png` | 333851 | 280483 |
| `README_supplementary.txt` | 1554 | 738 |

Total ZIP size: **743747 bytes (0.709 MiB)**. This is below the requested 250 MB limit.

## Uncertainty / limitations

- Table S4's confidence intervals are the already-stored 1,000-resample percentile intervals from the historical single test split; they were not recalculated.
- Table S5 combines the already-stored five-fold Wilcoxon p-values with the already-stored paired single-test-split bootstrap differences. The source explicitly notes low Wilcoxon power at `n = 5`.
- No missing requested source value or caption was encountered.
