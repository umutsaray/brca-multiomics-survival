# Phase 10 exploratory fold-stratified comparisons and DeepSurv code audit

Method: Harrell C within each estimable outer fold; delta = configuration minus Clinical-only; unweighted mean across estimable folds; participant bootstrap independently within folds; 9,999 iterations; seed 42; percentile 95% CI; two-sided empirical tail p-value with +1 finite-sample correction; BH-FDR across 24 tests separately per cohort.

## Old pooled versus Phase 10 fold-stratified results

| Cohort | Model | Configuration | Folds | Old pooled Δ | Old 95% CI | Old p | Old q | New fold-mean Δ | New 95% CI | New p | New q |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| mixed | coxph | RNA-seq only | 1,2,3,4,5 (5/5) | -0.156396 | [-0.243341, -0.071256] | 0.001998 | 0.011988 | -0.161264 | [-0.248549, -0.071853] | 0.000800 | 0.003200 |
| mixed | coxph | Methylation only | 1,2,3,4,5 (5/5) | -0.132681 | [-0.218686, -0.051664] | 0.001998 | 0.011988 | -0.135024 | [-0.221018, -0.049617] | 0.003600 | 0.010800 |
| mixed | coxph | miRNA only | 1,2,3,4,5 (5/5) | -0.179165 | [-0.266992, -0.091217] | 0.001998 | 0.011988 | -0.165865 | [-0.258616, -0.073596] | 0.000200 | 0.002400 |
| mixed | coxph | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.072949 | [-0.144230, -0.003916] | 0.039960 | 0.106560 | -0.070301 | [-0.142910, 0.001425] | 0.054800 | 0.082200 |
| mixed | coxph | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.071014 | [-0.143720, 0.000472] | 0.055944 | 0.134266 | -0.066988 | [-0.141541, 0.006043] | 0.074200 | 0.101600 |
| mixed | coxph | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.078426 | [-0.147677, -0.011571] | 0.021978 | 0.071928 | -0.079149 | [-0.155154, -0.006143] | 0.033000 | 0.052800 |
| mixed | coxph | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.035431 | [-0.102657, 0.027061] | 0.255744 | 0.306893 | -0.041020 | [-0.104734, 0.019918] | 0.197400 | 0.205983 |
| mixed | coxph | Clinical + miRNA | 1,2,3,4,5 (5/5) | -0.030997 | [-0.081240, 0.018061] | 0.195804 | 0.266400 | -0.027730 | [-0.081012, 0.023805] | 0.283200 | 0.283200 |
| mixed | rsf | RNA-seq only | 1,2,3,4,5 (5/5) | -0.119411 | [-0.204172, -0.037427] | 0.007992 | 0.038362 | -0.127988 | [-0.212219, -0.037572] | 0.005000 | 0.013333 |
| mixed | rsf | Methylation only | 1,2,3,4,5 (5/5) | -0.097000 | [-0.168605, -0.023875] | 0.011988 | 0.047952 | -0.108645 | [-0.188831, -0.021918] | 0.015600 | 0.028800 |
| mixed | rsf | miRNA only | 1,2,3,4,5 (5/5) | -0.164971 | [-0.254380, -0.074142] | 0.001998 | 0.011988 | -0.167211 | [-0.253043, -0.073670] | 0.000600 | 0.002880 |
| mixed | rsf | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.055222 | [-0.124840, 0.011822] | 0.103896 | 0.195497 | -0.060054 | [-0.127436, 0.011256] | 0.098600 | 0.118320 |
| mixed | rsf | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.057331 | [-0.124828, 0.013171] | 0.105894 | 0.195497 | -0.068067 | [-0.142861, 0.011290] | 0.091400 | 0.115453 |
| mixed | rsf | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.052375 | [-0.121854, 0.016243] | 0.137862 | 0.236335 | -0.057925 | [-0.129863, 0.017833] | 0.134000 | 0.153143 |
| mixed | rsf | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.060787 | [-0.124410, 0.006384] | 0.081918 | 0.178730 | -0.069633 | [-0.142051, 0.007971] | 0.076200 | 0.101600 |
| mixed | rsf | Clinical + miRNA | 1,2,3,4,5 (5/5) | -0.039898 | [-0.101860, 0.020583] | 0.189810 | 0.266400 | -0.044887 | [-0.109009, 0.023225] | 0.185400 | 0.202255 |
| mixed | deepsurv | RNA-seq only | 1,2,3,4,5 (5/5) | -0.033931 | [-0.141601, 0.072883] | 0.507493 | 0.553628 | -0.136962 | [-0.234032, -0.040094] | 0.006200 | 0.013527 |
| mixed | deepsurv | Methylation only | 3,4,5 (3/5) | -0.032381 | [-0.134509, 0.061829] | 0.507493 | 0.553628 | -0.188150 | [-0.303857, -0.066019] | 0.005600 | 0.013440 |
| mixed | deepsurv | miRNA only | 1,2,3,4,5 (5/5) | -0.122682 | [-0.229576, -0.025347] | 0.023976 | 0.071928 | -0.215130 | [-0.315423, -0.114781] | 0.000200 | 0.002400 |
| mixed | deepsurv | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.058515 | [-0.150166, 0.025293] | 0.199800 | 0.266400 | -0.156021 | [-0.237977, -0.069853] | 0.000600 | 0.002880 |
| mixed | deepsurv | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.068036 | [-0.166040, 0.028403] | 0.189810 | 0.266400 | -0.170959 | [-0.255779, -0.082204] | 0.000400 | 0.002880 |
| mixed | deepsurv | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.012455 | [-0.106112, 0.076464] | 0.779221 | 0.813100 | -0.119378 | [-0.202469, -0.031513] | 0.009800 | 0.019600 |
| mixed | deepsurv | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.053701 | [-0.143681, 0.031885] | 0.251748 | 0.306893 | -0.111380 | [-0.184779, -0.038479] | 0.002400 | 0.008229 |
| mixed | deepsurv | Clinical + miRNA | 1,2,3,4,5 (5/5) | 0.006912 | [-0.072225, 0.089068] | 0.891109 | 0.891109 | -0.073487 | [-0.139260, -0.007580] | 0.030600 | 0.052457 |
| primary | coxph | RNA-seq only | 1,2,3,4,5 (5/5) | -0.208899 | [-0.292934, -0.118682] | 0.001998 | 0.003996 | -0.204673 | [-0.286885, -0.111531] | 0.000200 | 0.000480 |
| primary | coxph | Methylation only | 2,4,5 (3/5) | -0.149365 | [-0.264727, -0.030652] | 0.017982 | 0.022714 | -0.165103 | [-0.287272, -0.040625] | 0.009600 | 0.014400 |
| primary | coxph | miRNA only | 1,2,3,4,5 (5/5) | -0.268078 | [-0.341862, -0.195268] | 0.001998 | 0.003996 | -0.259985 | [-0.325646, -0.186782] | 0.000200 | 0.000480 |
| primary | coxph | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.086109 | [-0.148183, -0.023891] | 0.017982 | 0.022714 | -0.083076 | [-0.149466, -0.012879] | 0.021400 | 0.027537 |
| primary | coxph | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.088478 | [-0.151342, -0.026255] | 0.009990 | 0.014985 | -0.088638 | [-0.153629, -0.019373] | 0.010800 | 0.015247 |
| primary | coxph | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.114469 | [-0.178747, -0.051938] | 0.001998 | 0.003996 | -0.101668 | [-0.165277, -0.035642] | 0.002400 | 0.003840 |
| primary | coxph | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.025637 | [-0.058392, 0.004883] | 0.127872 | 0.139497 | -0.024429 | [-0.065656, 0.014177] | 0.201800 | 0.220145 |
| primary | coxph | Clinical + miRNA | 1,2,3,4,5 (5/5) | -0.030359 | [-0.053954, -0.007144] | 0.009990 | 0.014985 | -0.028138 | [-0.053439, 0.001437] | 0.060200 | 0.068800 |
| primary | rsf | RNA-seq only | 1,2,3,4,5 (5/5) | -0.222386 | [-0.310644, -0.123361] | 0.001998 | 0.003996 | -0.203145 | [-0.288321, -0.111399] | 0.000200 | 0.000480 |
| primary | rsf | Methylation only | 2,4,5 (3/5) | -0.170663 | [-0.299350, -0.020501] | 0.019980 | 0.023976 | -0.183787 | [-0.335165, -0.021466] | 0.024000 | 0.028800 |
| primary | rsf | miRNA only | 1,2,3,4,5 (5/5) | -0.230337 | [-0.299117, -0.148872] | 0.001998 | 0.003996 | -0.229776 | [-0.295103, -0.158394] | 0.000200 | 0.000480 |
| primary | rsf | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.149211 | [-0.219107, -0.071048] | 0.001998 | 0.003996 | -0.135200 | [-0.211045, -0.054526] | 0.001200 | 0.002215 |
| primary | rsf | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.142044 | [-0.211959, -0.064209] | 0.001998 | 0.003996 | -0.126329 | [-0.200683, -0.049482] | 0.001400 | 0.002400 |
| primary | rsf | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.144074 | [-0.213223, -0.072581] | 0.001998 | 0.003996 | -0.127745 | [-0.199516, -0.052141] | 0.001200 | 0.002215 |
| primary | rsf | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.040110 | [-0.092980, 0.012848] | 0.141858 | 0.148026 | -0.038342 | [-0.100380, 0.021855] | 0.217000 | 0.226435 |
| primary | rsf | Clinical + miRNA | 1,2,3,4,5 (5/5) | -0.018286 | [-0.048582, 0.014257] | 0.285714 | 0.285714 | -0.009627 | [-0.046824, 0.027766] | 0.587600 | 0.587600 |
| primary | deepsurv | RNA-seq only | 1,2,3,4,5 (5/5) | -0.172788 | [-0.244014, -0.096515] | 0.001998 | 0.003996 | -0.211740 | [-0.299993, -0.116794] | 0.000200 | 0.000480 |
| primary | deepsurv | Methylation only | 1,5 (2/5) | -0.288597 | [-0.413680, -0.165460] | 0.001998 | 0.003996 | -0.346299 | [-0.470674, -0.203257] | 0.000200 | 0.000480 |
| primary | deepsurv | miRNA only | 1,2,3,4,5 (5/5) | -0.245532 | [-0.323178, -0.166095] | 0.001998 | 0.003996 | -0.243292 | [-0.314538, -0.165545] | 0.000200 | 0.000480 |
| primary | deepsurv | Clinical + RNA-seq | 1,2,3,4,5 (5/5) | -0.130494 | [-0.221154, -0.038912] | 0.011988 | 0.016924 | -0.147739 | [-0.220310, -0.064345] | 0.000600 | 0.001309 |
| primary | deepsurv | Clinical + RNA-seq + Methylation | 1,2,3,4,5 (5/5) | -0.124250 | [-0.199449, -0.041579] | 0.007992 | 0.013701 | -0.179419 | [-0.258134, -0.091253] | 0.000200 | 0.000480 |
| primary | deepsurv | Clinical + RNA-seq + miRNA | 1,2,3,4,5 (5/5) | -0.169773 | [-0.262507, -0.080458] | 0.003996 | 0.007377 | -0.216610 | [-0.299953, -0.124305] | 0.000200 | 0.000480 |
| primary | deepsurv | Clinical + Methylation | 1,2,3,4,5 (5/5) | -0.046231 | [-0.096519, 0.006666] | 0.079920 | 0.091337 | -0.057699 | [-0.103974, -0.008364] | 0.021800 | 0.027537 |
| primary | deepsurv | Clinical + miRNA | 1,2,3,4,5 (5/5) | -0.159500 | [-0.254934, -0.063264] | 0.001998 | 0.003996 | -0.164634 | [-0.238710, -0.082909] | 0.000200 | 0.000480 |

## Methylation-only estimable-fold warning

These rows use each model’s own estimable fold set. Consequently, methylation-only results across CoxPH, RSF, and DeepSurv are not directly comparable.

| Cohort | Model | Number of folds | Fold IDs |
|---|---|---:|---|
| mixed | coxph | 5 | 1, 2, 3, 4, 5 |
| mixed | rsf | 5 | 1, 2, 3, 4, 5 |
| mixed | deepsurv | 3 | 3, 4, 5 |
| primary | coxph | 3 | 2, 4, 5 |
| primary | rsf | 3 | 2, 4, 5 |
| primary | deepsurv | 2 | 1, 5 |

## Rows with direction or significance changes

- mixed / coxph / Clinical + RNA-seq: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=no; old Δ=-0.072949, p=0.039960, q=0.106560; new Δ=-0.070301, p=0.054800, q=0.082200.
- mixed / deepsurv / RNA-seq only: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.033931, p=0.507493, q=0.553628; new Δ=-0.136962, p=0.006200, q=0.013527.
- mixed / deepsurv / Methylation only: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.032381, p=0.507493, q=0.553628; new Δ=-0.188150, p=0.005600, q=0.013440.
- mixed / deepsurv / miRNA only: direction changed=no; raw-CI significance changed=no; BH q<0.05 status changed=yes; old Δ=-0.122682, p=0.023976, q=0.071928; new Δ=-0.215130, p=0.000200, q=0.002400.
- mixed / deepsurv / Clinical + RNA-seq: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.058515, p=0.199800, q=0.266400; new Δ=-0.156021, p=0.000600, q=0.002880.
- mixed / deepsurv / Clinical + RNA-seq + Methylation: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.068036, p=0.189810, q=0.266400; new Δ=-0.170959, p=0.000400, q=0.002880.
- mixed / deepsurv / Clinical + RNA-seq + miRNA: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.012455, p=0.779221, q=0.813100; new Δ=-0.119378, p=0.009800, q=0.019600.
- mixed / deepsurv / Clinical + Methylation: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.053701, p=0.251748, q=0.306893; new Δ=-0.111380, p=0.002400, q=0.008229.
- mixed / deepsurv / Clinical + miRNA: direction changed=yes; raw-CI significance changed=yes; BH q<0.05 status changed=no; old Δ=0.006912, p=0.891109, q=0.891109; new Δ=-0.073487, p=0.030600, q=0.052457.
- primary / coxph / Clinical + miRNA: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.030359, p=0.009990, q=0.014985; new Δ=-0.028138, p=0.060200, q=0.068800.
- primary / deepsurv / Clinical + Methylation: direction changed=no; raw-CI significance changed=yes; BH q<0.05 status changed=yes; old Δ=-0.046231, p=0.079920, q=0.091337; new Δ=-0.057699, p=0.021800, q=0.027537.

## Code-verified DeepSurv procedure

The premise that the network is retrained after selecting an epoch is not implemented in the retained Phase 7/9 code. Exact answers:

1. **No reinitialization/retraining occurs.** `fit_network` trains one network while monitoring internal-validation C-index, stores a copy of the best state dictionary, and after stopping loads that checkpoint (`src/models/phase7_deepsurv.py`, `fit_network`, lines 31–40).
2. **Optimization uses internal-training only.** The loss tensors and feature blocks used by the optimizer are constructed from `xf`/`yf`. Internal validation `xv`/`yv` is used to compute the early-stopping C-index; it is not merged back into training (lines 31–40).
3. **No preprocessing refit follows epoch selection.** Phase 7 fits clinical encoding and methylation imputation on internal train (lines 21–29), and performs internal-train outcome-guided selection before the fit (lines 50–57). Validation/test are transform-only. Phase 9 loads those cached internal-train-derived frames and passes them to the same fit function (`src/models/phase9_deepsurv.py`, lines 31–42).
4. **Outer-test risk comes from the restored best checkpoint.** After `model.load_state_dict(best)`, `predict_risk(model, tb)` produces the outer-test risk. It is not produced by a separately retrained model (`src/models/phase7_deepsurv.py`, line 40).

Therefore manuscript wording that says the network “is retrained” after epoch selection is not code-accurate. The code performs checkpoint selection and restoration only.

## Seeds

- Outer folds: seed 42.
- Internal participant-grouped/event-stratified split: 4200 + outer fold.
- Initial and only DeepSurv training run: 42 + outer fold.
- Phase 10 fold-stratified bootstrap: seed 42, 9,999 iterations.
