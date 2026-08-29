# Phase 7 updated additional metrics

Only DeepSurv OOF risks were replaced. CoxPH and RSF were retained and verified exactly against the preceding results.

Maximum absolute CoxPH/RSF metric difference: `0.0`.

## mixed

| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |
|---|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | coxph | 0.6279 | nan | 0.7248 | 0.7107 | nan | 0.0888 |
| Clinical-only | rsf | 0.6295 | nan | 0.7328 | 0.6907 | nan | 0.0803 |
| Clinical-only | deepsurv | 0.6018 | nan | 0.7307 | 0.7168 | nan | 0.0894 |
| Full | coxph | 0.6472 | nan | 0.7173 | 0.6283 | nan | 0.0922 |
| Full | rsf | 0.6563 | nan | 0.6887 | 0.6592 | nan | 0.0889 |
| Full | deepsurv | 0.6514 | nan | 0.6332 | 0.5944 | nan | 0.0913 |

## primary

| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |
|---|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | coxph | 0.5310 | 0.8914 | 0.7769 | 0.7322 | 0.7614 | 0.0813 |
| Clinical-only | rsf | 0.5325 | 0.8839 | 0.7666 | 0.7168 | 0.7497 | 0.0736 |
| Clinical-only | deepsurv | 0.5210 | 0.8860 | 0.7611 | 0.7065 | 0.7388 | 0.0810 |
| Full | coxph | 0.5147 | 0.7342 | 0.6623 | 0.6278 | 0.6511 | 0.0833 |
| Full | rsf | 0.5135 | 0.6490 | 0.6650 | 0.6291 | 0.6254 | 0.0827 |
| Full | deepsurv | 0.4461 | 0.4856 | 0.5031 | 0.5771 | 0.5470 | 0.0835 |
