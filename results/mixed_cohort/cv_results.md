# Five-Fold Stratified Cross-Validation Results

Samples: 842; folds: 5; stratification: event; seed: 42.
Feature selection was recomputed using only each outer fold's training subset.

| Fold | Train/Test | CoxPH C-index | RSF C-index | DeepSurv C-index | Selection s | Total s |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 673/169 | 0.656944 | 0.655675 | 0.621433 | 75.80 | 173.35 |
| 2 | 673/169 | 0.647088 | 0.787266 | 0.675222 | 75.40 | 98.78 |
| 3 | 674/168 | 0.741935 | 0.646628 | 0.624633 | 75.57 | 106.34 |
| 4 | 674/168 | 0.593829 | 0.509446 | 0.515743 | 72.75 | 96.21 |
| 5 | 674/168 | 0.588426 | 0.683613 | 0.619794 | 75.67 | 106.09 |

## Mean ± sample standard deviation

- CoxPH: 0.645644 ± 0.061969
- Random Survival Forest: 0.656525 ± 0.099446
- DeepSurv: 0.611365 ± 0.058245

Total elapsed time: 580.78 seconds.
