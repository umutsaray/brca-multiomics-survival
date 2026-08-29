# Fold-stratified clinical-only versus full multi-omics comparison

The new analysis resamples participants independently within every outer fold, never compares raw risks across folds, computes Harrell C within each fold, and uses the unweighted mean of five fold-specific deltas. Delta is Clinical-only minus Full; positive values favor Clinical-only. All models share the same participant resample indices in each of 9,999 iterations (seed 42).

## Old pooled versus new fold-stratified results

The old stored estimates used Full minus Clinical. They are sign-normalized below to Clinical minus Full for a direct comparison.

| Cohort | Model | Old pooled delta | Old 95% CI | Old p | New fold-stratified delta | New 95% CI | New p |
|---|---|---:|---:|---:|---:|---:|---:|
| mixed | coxph | 0.0682 | [0.0019, 0.1373] | 0.0448 | 0.0655 | [-0.0049, 0.1352] | 0.0642 |
| mixed | rsf | 0.0609 | [-0.0110, 0.1309] | 0.095 | 0.0711 | [-0.0093, 0.1466] | 0.0816 |
| mixed | deepsurv | 0.0161 | [-0.0836, 0.1200] | 0.7692 | 0.1474 | [0.0518, 0.2388] | 0.0034 |
| primary | coxph | 0.1229 | [0.0550, 0.1884] | 0.0002 | 0.1073 | [0.0370, 0.1727] | 0.0034 |
| primary | rsf | 0.1473 | [0.0735, 0.2166] | 0.0004 | 0.1352 | [0.0584, 0.2092] | 0.0008 |
| primary | deepsurv | 0.2149 | [0.1518, 0.2776] | 0.0002 | 0.2180 | [0.1194, 0.3051] | 0.0002 |

## Fold-level new results

### mixed

| Model | Fold | Clinical C | Full C | Delta |
|---|---:|---:|---:|---:|
| coxph | 1 | 0.6677 | 0.6210 | 0.0467 |
| coxph | 2 | 0.7495 | 0.7548 | -0.0054 |
| coxph | 3 | 0.7145 | 0.6601 | 0.0543 |
| coxph | 4 | 0.7923 | 0.7146 | 0.0777 |
| coxph | 5 | 0.7549 | 0.6007 | 0.1542 |
| rsf | 1 | 0.6673 | 0.5588 | 0.1086 |
| rsf | 2 | 0.6984 | 0.7806 | -0.0822 |
| rsf | 3 | 0.7841 | 0.6087 | 0.1754 |
| rsf | 4 | 0.7850 | 0.7163 | 0.0687 |
| rsf | 5 | 0.7589 | 0.6738 | 0.0852 |
| deepsurv | 1 | 0.6833 | 0.5790 | 0.1043 |
| deepsurv | 2 | 0.7052 | 0.6994 | 0.0058 |
| deepsurv | 3 | 0.7362 | 0.4326 | 0.3036 |
| deepsurv | 4 | 0.7811 | 0.6303 | 0.1508 |
| deepsurv | 5 | 0.7917 | 0.6191 | 0.1726 |

### primary

| Model | Fold | Clinical C | Full C | Delta |
|---|---:|---:|---:|---:|
| coxph | 1 | 0.8271 | 0.7402 | 0.0869 |
| coxph | 2 | 0.6363 | 0.6615 | -0.0253 |
| coxph | 3 | 0.7958 | 0.6200 | 0.1758 |
| coxph | 4 | 0.7361 | 0.6814 | 0.0547 |
| coxph | 5 | 0.8037 | 0.5592 | 0.2446 |
| rsf | 1 | 0.8041 | 0.6428 | 0.1614 |
| rsf | 2 | 0.6533 | 0.6538 | -0.0005 |
| rsf | 3 | 0.8178 | 0.5987 | 0.2191 |
| rsf | 4 | 0.7348 | 0.6565 | 0.0783 |
| rsf | 5 | 0.7762 | 0.5583 | 0.2179 |
| deepsurv | 1 | 0.8298 | 0.5690 | 0.2608 |
| deepsurv | 2 | 0.5890 | 0.4857 | 0.1033 |
| deepsurv | 3 | 0.7852 | 0.5426 | 0.2425 |
| deepsurv | 4 | 0.7064 | 0.5457 | 0.1607 |
| deepsurv | 5 | 0.8113 | 0.4883 | 0.3229 |

## Direction and significance

- mixed, coxph: delta 0.0655, 95% CI [-0.0049, 0.1352], p=0.0642; the CI includes zero; no statistically detectable difference.
- mixed, rsf: delta 0.0711, 95% CI [-0.0093, 0.1466], p=0.0816; the CI includes zero; no statistically detectable difference.
- mixed, deepsurv: delta 0.1474, 95% CI [0.0518, 0.2388], p=0.0034; Clinical-only is favored and the CI excludes zero.
- primary, coxph: delta 0.1073, 95% CI [0.0370, 0.1727], p=0.0034; Clinical-only is favored and the CI excludes zero.
- primary, rsf: delta 0.1352, 95% CI [0.0584, 0.2092], p=0.0008; Clinical-only is favored and the CI excludes zero.
- primary, deepsurv: delta 0.2180, 95% CI [0.1194, 0.3051], p=0.0002; Clinical-only is favored and the CI excludes zero.

The prior pooled values are retained only for audit comparison; they are not the Phase 7 primary estimates.