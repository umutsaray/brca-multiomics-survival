# Phase 9 remaining DeepSurv ablation correction

Requested fit attempts: 80; completed: 75; skipped because internal-train FDR selection yielded zero requested-block features: 5.

No fallback features and no old/leaky predictions were substituted for skipped folds. Thus every retained DeepSurv risk in the current ten-configuration grid comes from the leakage-corrected pipeline.

## Old versus corrected DeepSurv C-index

### mixed

| Configuration | Fold | Old C | Corrected C | Corrected − old |
|---|---:|---:|---:|---:|
| RNA-seq only | 1 | 0.454475 | 0.613230 | 0.158755 |
| RNA-seq only | 2 | 0.660369 | 0.667239 | 0.006870 |
| RNA-seq only | 3 | 0.592754 | 0.573913 | -0.018841 |
| RNA-seq only | 4 | 0.581087 | 0.573813 | -0.007274 |
| RNA-seq only | 5 | 0.647871 | 0.584580 | -0.063291 |
| Methylation only | 1 | 0.527626 | NA | NA |
| Methylation only | 2 | 0.654787 | NA | NA |
| Methylation only | 3 | 0.563043 | 0.584783 | 0.021739 |
| Methylation only | 4 | 0.596063 | 0.623021 | 0.026958 |
| Methylation only | 5 | 0.562716 | 0.536824 | -0.025892 |
| miRNA only | 1 | 0.561868 | 0.457588 | -0.104280 |
| miRNA only | 2 | 0.639330 | 0.586088 | -0.053242 |
| miRNA only | 3 | 0.501449 | 0.376087 | -0.125362 |
| miRNA only | 4 | 0.561831 | 0.516902 | -0.044929 |
| miRNA only | 5 | 0.478999 | 0.685270 | 0.206272 |
| Clinical + RNA-seq | 1 | 0.431128 | 0.428794 | -0.002335 |
| Clinical + RNA-seq | 2 | 0.690854 | 0.660799 | -0.030056 |
| Clinical + RNA-seq | 3 | 0.621014 | 0.576087 | -0.044928 |
| Clinical + RNA-seq | 4 | 0.615747 | 0.618314 | 0.002567 |
| Clinical + RNA-seq | 5 | 0.668009 | 0.633487 | -0.034522 |
| Clinical + RNA-seq + Methylation | 1 | 0.511284 | 0.428794 | -0.082490 |
| Clinical + RNA-seq + Methylation | 2 | 0.701589 | 0.660799 | -0.040790 |
| Clinical + RNA-seq + Methylation | 3 | 0.643478 | 0.526812 | -0.116667 |
| Clinical + RNA-seq + Methylation | 4 | 0.606761 | 0.630295 | 0.023534 |
| Clinical + RNA-seq + Methylation | 5 | 0.669160 | 0.596087 | -0.073072 |
| Clinical + RNA-seq + miRNA | 1 | 0.560311 | 0.578988 | 0.018677 |
| Clinical + RNA-seq + miRNA | 2 | 0.686990 | 0.699442 | 0.012452 |
| Clinical + RNA-seq + miRNA | 3 | 0.640580 | 0.555072 | -0.085507 |
| Clinical + RNA-seq + miRNA | 4 | 0.711596 | 0.618742 | -0.092854 |
| Clinical + RNA-seq + miRNA | 5 | 0.678366 | 0.648446 | -0.029919 |
| Clinical + Methylation | 1 | 0.678599 | 0.683268 | 0.004669 |
| Clinical + Methylation | 2 | 0.692572 | 0.705238 | 0.012666 |
| Clinical + Methylation | 3 | 0.499275 | 0.632609 | 0.133333 |
| Clinical + Methylation | 4 | 0.606761 | 0.606333 | -0.000428 |
| Clinical + Methylation | 5 | 0.628884 | 0.513234 | -0.115650 |
| Clinical + miRNA | 1 | 0.593774 | 0.600778 | 0.007004 |
| Clinical + miRNA | 2 | 0.642336 | 0.555174 | -0.087162 |
| Clinical + miRNA | 3 | 0.666667 | 0.673188 | 0.006522 |
| Clinical + miRNA | 4 | 0.580231 | 0.796748 | 0.216517 |
| Clinical + miRNA | 5 | 0.770426 | 0.704258 | -0.066168 |

Mean ± sample SD (means use only folds with an estimable model):

| Configuration | Old | Corrected |
|---|---:|---:|
| RNA-seq only | 0.587311 ± 0.081730 (n=5 folds) | 0.602555 ± 0.039584 (n=5 folds) |
| Methylation only | 0.580847 ± 0.047898 (n=5 folds) | 0.581543 ± 0.043190 (n=3 folds) |
| miRNA only | 0.548695 ± 0.062551 (n=5 folds) | 0.524387 ± 0.118534 (n=5 folds) |
| Clinical + RNA-seq | 0.605351 ± 0.102401 (n=5 folds) | 0.583496 ± 0.091753 (n=5 folds) |
| Clinical + RNA-seq + Methylation | 0.626454 ± 0.073162 (n=5 folds) | 0.568557 ± 0.092675 (n=5 folds) |
| Clinical + RNA-seq + miRNA | 0.655569 ± 0.059044 (n=5 folds) | 0.620138 ± 0.057051 (n=5 folds) |
| Clinical + Methylation | 0.621218 ± 0.076684 (n=5 folds) | 0.628136 ± 0.075303 (n=5 folds) |
| Clinical + miRNA | 0.650687 ± 0.075607 (n=5 folds) | 0.666029 ± 0.093737 (n=5 folds) |

### primary

| Configuration | Fold | Old C | Corrected C | Corrected − old |
|---|---:|---:|---:|---:|
| RNA-seq only | 1 | 0.520298 | 0.460758 | -0.059540 |
| RNA-seq only | 2 | 0.647253 | 0.458242 | -0.189011 |
| RNA-seq only | 3 | 0.422585 | 0.495028 | 0.072443 |
| RNA-seq only | 4 | 0.472299 | 0.546399 | 0.074100 |
| RNA-seq only | 5 | 0.441667 | 0.702500 | 0.260833 |
| Methylation only | 1 | NA | 0.499323 | NA |
| Methylation only | 2 | 0.521978 | NA | NA |
| Methylation only | 3 | NA | NA | NA |
| Methylation only | 4 | 0.639197 | NA | NA |
| Methylation only | 5 | 0.530833 | 0.449167 | -0.081667 |
| miRNA only | 1 | 0.500000 | 0.604195 | 0.104195 |
| miRNA only | 2 | 0.428571 | 0.454396 | 0.025824 |
| miRNA only | 3 | 0.492898 | 0.443537 | -0.049361 |
| miRNA only | 4 | 0.500000 | 0.505540 | 0.005540 |
| miRNA only | 5 | 0.544583 | 0.497500 | -0.047083 |
| Clinical + RNA-seq | 1 | 0.548038 | 0.560217 | 0.012179 |
| Clinical + RNA-seq | 2 | 0.707692 | 0.507692 | -0.200000 |
| Clinical + RNA-seq | 3 | 0.563210 | 0.578125 | 0.014915 |
| Clinical + RNA-seq | 4 | 0.479917 | 0.617729 | 0.137812 |
| Clinical + RNA-seq | 5 | 0.532917 | 0.719167 | 0.186250 |
| Clinical + RNA-seq + Methylation | 1 | 0.548038 | 0.578484 | 0.030447 |
| Clinical + RNA-seq + Methylation | 2 | 0.620879 | 0.507692 | -0.113187 |
| Clinical + RNA-seq + Methylation | 3 | 0.563210 | 0.578125 | 0.014915 |
| Clinical + RNA-seq + Methylation | 4 | 0.590028 | 0.617729 | 0.027701 |
| Clinical + RNA-seq + Methylation | 5 | 0.502500 | 0.542500 | 0.040000 |
| Clinical + RNA-seq + miRNA | 1 | 0.594723 | 0.559540 | -0.035183 |
| Clinical + RNA-seq + miRNA | 2 | 0.536264 | 0.485714 | -0.050549 |
| Clinical + RNA-seq + miRNA | 3 | 0.470170 | 0.542614 | 0.072443 |
| Clinical + RNA-seq + miRNA | 4 | 0.448753 | 0.545706 | 0.096953 |
| Clinical + RNA-seq + miRNA | 5 | 0.520833 | 0.505000 | -0.015833 |
| Clinical + Methylation | 1 | 0.824763 | 0.801759 | -0.023004 |
| Clinical + Methylation | 2 | 0.580220 | 0.589011 | 0.008791 |
| Clinical + Methylation | 3 | 0.816051 | 0.785156 | -0.030895 |
| Clinical + Methylation | 4 | 0.648892 | 0.706371 | 0.057479 |
| Clinical + Methylation | 5 | 0.505833 | 0.550833 | 0.045000 |
| Clinical + miRNA | 1 | 0.837957 | 0.510149 | -0.327808 |
| Clinical + miRNA | 2 | 0.569231 | 0.553297 | -0.015934 |
| Clinical + miRNA | 3 | 0.758523 | 0.515980 | -0.242543 |
| Clinical + miRNA | 4 | 0.569945 | 0.549861 | -0.020083 |
| Clinical + miRNA | 5 | 0.767917 | 0.769167 | 0.001250 |

Mean ± sample SD (means use only folds with an estimable model):

| Configuration | Old | Corrected |
|---|---:|---:|
| RNA-seq only | 0.500820 ± 0.089799 (n=5 folds) | 0.532585 ± 0.101438 (n=5 folds) |
| Methylation only | 0.564003 ± 0.065270 (n=3 folds) | 0.474245 ± 0.035466 (n=2 folds) |
| miRNA only | 0.493210 ± 0.041562 (n=5 folds) | 0.501034 ± 0.063553 (n=5 folds) |
| Clinical + RNA-seq | 0.566355 ± 0.085017 (n=5 folds) | 0.596586 ± 0.079118 (n=5 folds) |
| Clinical + RNA-seq + Methylation | 0.564931 ± 0.044570 (n=5 folds) | 0.564906 ± 0.041608 (n=5 folds) |
| Clinical + RNA-seq + miRNA | 0.514149 ± 0.057525 (n=5 folds) | 0.527715 ± 0.030978 (n=5 folds) |
| Clinical + Methylation | 0.675152 ± 0.141957 (n=5 folds) | 0.686626 ± 0.113272 (n=5 folds) |
| Clinical + miRNA | 0.700714 ± 0.123574 (n=5 folds) | 0.579691 ± 0.107683 (n=5 folds) |

## Direction and significance audit

Table 4/4b are descriptive and contain no inferential significance test; significance change is therefore not defined for those tables. Direction below refers to the pooled Supplementary Table S2 delta (configuration minus corrected Clinical-only). BH significance uses q < 0.05.

| Cohort | Configuration | Old delta | New delta | Direction changed? | Raw-CI significance changed? | BH significance changed? |
|---|---|---:|---:|---|---|---|
| mixed | RNA-seq only | -0.039594 | -0.033931 | No | No | No |
| mixed | Methylation only | -0.048506 | -0.032381 | No | No | No |
| mixed | miRNA only | -0.044930 | -0.122682 | No | Yes | No |
| mixed | Clinical + RNA-seq | 0.036268 | -0.058515 | Yes | No | No |
| mixed | Clinical + RNA-seq + Methylation | 0.035051 | -0.068036 | Yes | No | No |
| mixed | Clinical + RNA-seq + miRNA | 0.047071 | -0.012455 | Yes | No | No |
| mixed | Clinical + Methylation | -0.041289 | -0.053701 | No | No | No |
| mixed | Clinical + miRNA | 0.044376 | 0.006912 | No | No | No |
| primary | RNA-seq only | -0.207561 | -0.172788 | No | No | No |
| primary | Methylation only | -0.179502 | -0.288597 | No | No | No |
| primary | miRNA only | -0.262911 | -0.245532 | No | No | No |
| primary | Clinical + RNA-seq | -0.179801 | -0.130494 | No | No | No |
| primary | Clinical + RNA-seq + Methylation | -0.210206 | -0.124250 | No | No | No |
| primary | Clinical + RNA-seq + miRNA | -0.245486 | -0.169773 | No | No | No |
| primary | Clinical + Methylation | -0.070776 | -0.046231 | No | Yes | Yes |
| primary | Clinical + miRNA | -0.090093 | -0.159500 | No | No | No |

## SHAP / Table 10 impact check

Table 10, Supplementary Table S1, and Figure S2 use only the primary-tumor-only Full DeepSurv models. Full was already retrained under the corrected pipeline in Phase 7 and is not one of the eight Phase 9 configurations. No SHAP recomputation is required; the Phase 7 SHAP outputs remain current.

## Supplementary Table S2 impact

The retained file contains 48 exploratory rows: 24 per cohort (8 configurations × 3 models). Exactly 16 rows have changed model inputs: the 8 DeepSurv exploratory comparisons in each of two cohorts. CoxPH and RSF inputs are unchanged, but their BH-adjusted q-values were recomputed jointly because the family includes the changed DeepSurv p-values. The complete updated table is `supplementary/table_s2_exploratory_clinical_vs_omics.csv`; full raw arrays are in `clinical_vs_omics_phase9.json`.

## Method and seeds

- Outer split: five-fold participant-grouped/event-stratified StratifiedGroupKFold, seed 42.
- Internal split: first held fold of seven-fold participant-grouped/event-stratified StratifiedGroupKFold, seed 4200 + outer fold.
- Feature selection, imputation, clinical encoding, and scaling: fitted on internal train only.
- DeepSurv training seed: 42 + outer fold; validation used only for early stopping; outer test predicted once.
- Network/optimizer/loss settings are identical to Phase 7.
- Supplementary Table S2: pooled participant bootstrap, 1,000 iterations, seed 42, BH correction over 24 exploratory comparisons within each cohort.

Attempted 80 new fits; completed 75; skipped 5 zero-feature Methylation-only folds. Together with the 20 Phase 7 fits, 100 corrected fit attempts cover all ten configurations; 95 models were estimable and fitted.
Generation date: 2026-08-29. Git provenance: the commit containing this report; its immutable hash is reported in the delivery message because a commit cannot embed its own hash.