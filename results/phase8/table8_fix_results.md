# Table 8 Phase 8 correction

Mixed cohort, Full multi-omics configuration. DeepSurv risks are the Phase 7 leakage-corrected out-of-fold predictions; CoxPH and RSF risks are the unchanged retained Full-configuration predictions.

New method: Harrell C is calculated separately inside each outer fold. Fold delta is DeepSurv C minus comparator C. The point estimate is the unweighted mean of five fold deltas. In each of 9,999 iterations, participants are resampled with replacement independently inside every fold; all three models use identical resample indices. Two-sided p-values use the paired bootstrap distribution with a +1 finite-sample correction. Seed: 42.

## Old pooled versus current fold-stratified results

| Comparison | Old pooled/leaky delta | Old 95% CI | Old p | New fold-stratified/corrected delta | Fold SD | New 95% CI | New p | Direction changed? | Significance changed? |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| DeepSurv − COXPH | 0.002087 | [-0.063285, 0.061672] | 0.919081 | -0.078166 | 0.091517 | [-0.148998, -0.001643] | 0.046800 | Yes | Yes |
| DeepSurv − RSF | -0.035605 | [-0.093576, 0.017061] | 0.185814 | -0.075534 | 0.070461 | [-0.138978, -0.008169] | 0.028200 | No | Yes |

## Fold-level current values

| Fold | DeepSurv C | CoxPH C | RSF C | DeepSurv − CoxPH | DeepSurv − RSF |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.578988 | 0.621012 | 0.558755 | -0.042023 | 0.020233 |
| 2 | 0.699442 | 0.754830 | 0.780593 | -0.055389 | -0.081151 |
| 3 | 0.432609 | 0.660145 | 0.608696 | -0.227536 | -0.176087 |
| 4 | 0.630295 | 0.714591 | 0.716303 | -0.084296 | -0.086008 |
| 5 | 0.619102 | 0.600690 | 0.673763 | 0.018412 | -0.054661 |

## Observed change

- DeepSurv − COXPH: direction changed; significance status changed. Old and new significance statuses are not significant and significant, respectively.
- DeepSurv − RSF: direction did not change; significance status changed. Old and new significance statuses are not significant and significant, respectively.