# Phase 9 clinical-only versus omics pooled participant bootstrap

This intentionally retains the prior pooled-score method for Supplementary Table S2. Inputs are Phase 9 leakage-corrected DeepSurv risks for all retained configurations and unchanged CoxPH/RSF risks. Iterations: 1,000; seed: 42.

## mixed

| Role | Model | Configuration | Delta (configuration − Clinical) | 95% CI | p | BH-FDR p |
|---|---|---|---:|---:|---:|---:|
| exploratory | coxph | RNA-seq only | -0.156396 | [-0.243341, -0.071256] | 0.001998 | 0.011988 |
| exploratory | coxph | Methylation only | -0.132681 | [-0.218686, -0.051664] | 0.001998 | 0.011988 |
| exploratory | coxph | miRNA only | -0.179165 | [-0.266992, -0.091217] | 0.001998 | 0.011988 |
| exploratory | coxph | Clinical + RNA-seq | -0.072949 | [-0.144230, -0.003916] | 0.039960 | 0.106560 |
| exploratory | coxph | Clinical + RNA-seq + Methylation | -0.071014 | [-0.143720, 0.000472] | 0.055944 | 0.134266 |
| exploratory | coxph | Clinical + RNA-seq + miRNA | -0.078426 | [-0.147677, -0.011571] | 0.021978 | 0.071928 |
| exploratory | coxph | Clinical + Methylation | -0.035431 | [-0.102657, 0.027061] | 0.255744 | 0.306893 |
| exploratory | coxph | Clinical + miRNA | -0.030997 | [-0.081240, 0.018061] | 0.195804 | 0.266400 |
| exploratory | rsf | RNA-seq only | -0.119411 | [-0.204172, -0.037427] | 0.007992 | 0.038362 |
| exploratory | rsf | Methylation only | -0.097000 | [-0.168605, -0.023875] | 0.011988 | 0.047952 |
| exploratory | rsf | miRNA only | -0.164971 | [-0.254380, -0.074142] | 0.001998 | 0.011988 |
| exploratory | rsf | Clinical + RNA-seq | -0.055222 | [-0.124840, 0.011822] | 0.103896 | 0.195497 |
| exploratory | rsf | Clinical + RNA-seq + Methylation | -0.057331 | [-0.124828, 0.013171] | 0.105894 | 0.195497 |
| exploratory | rsf | Clinical + RNA-seq + miRNA | -0.052375 | [-0.121854, 0.016243] | 0.137862 | 0.236335 |
| exploratory | rsf | Clinical + Methylation | -0.060787 | [-0.124410, 0.006384] | 0.081918 | 0.178730 |
| exploratory | rsf | Clinical + miRNA | -0.039898 | [-0.101860, 0.020583] | 0.189810 | 0.266400 |
| exploratory | deepsurv | RNA-seq only | -0.033931 | [-0.141601, 0.072883] | 0.507493 | 0.553628 |
| exploratory | deepsurv | Methylation only | -0.032381 | [-0.134509, 0.061829] | 0.507493 | 0.553628 |
| exploratory | deepsurv | miRNA only | -0.122682 | [-0.229576, -0.025347] | 0.023976 | 0.071928 |
| exploratory | deepsurv | Clinical + RNA-seq | -0.058515 | [-0.150166, 0.025293] | 0.199800 | 0.266400 |
| exploratory | deepsurv | Clinical + RNA-seq + Methylation | -0.068036 | [-0.166040, 0.028403] | 0.189810 | 0.266400 |
| exploratory | deepsurv | Clinical + RNA-seq + miRNA | -0.012455 | [-0.106112, 0.076464] | 0.779221 | 0.813100 |
| exploratory | deepsurv | Clinical + Methylation | -0.053701 | [-0.143681, 0.031885] | 0.251748 | 0.306893 |
| exploratory | deepsurv | Clinical + miRNA | 0.006912 | [-0.072225, 0.089068] | 0.891109 | 0.891109 |

## primary

| Role | Model | Configuration | Delta (configuration − Clinical) | 95% CI | p | BH-FDR p |
|---|---|---|---:|---:|---:|---:|
| exploratory | coxph | RNA-seq only | -0.208899 | [-0.292934, -0.118682] | 0.001998 | 0.003996 |
| exploratory | coxph | Methylation only | -0.149365 | [-0.264727, -0.030652] | 0.017982 | 0.022714 |
| exploratory | coxph | miRNA only | -0.268078 | [-0.341862, -0.195268] | 0.001998 | 0.003996 |
| exploratory | coxph | Clinical + RNA-seq | -0.086109 | [-0.148183, -0.023891] | 0.017982 | 0.022714 |
| exploratory | coxph | Clinical + RNA-seq + Methylation | -0.088478 | [-0.151342, -0.026255] | 0.009990 | 0.014985 |
| exploratory | coxph | Clinical + RNA-seq + miRNA | -0.114469 | [-0.178747, -0.051938] | 0.001998 | 0.003996 |
| exploratory | coxph | Clinical + Methylation | -0.025637 | [-0.058392, 0.004883] | 0.127872 | 0.139497 |
| exploratory | coxph | Clinical + miRNA | -0.030359 | [-0.053954, -0.007144] | 0.009990 | 0.014985 |
| exploratory | rsf | RNA-seq only | -0.222386 | [-0.310644, -0.123361] | 0.001998 | 0.003996 |
| exploratory | rsf | Methylation only | -0.170663 | [-0.299350, -0.020501] | 0.019980 | 0.023976 |
| exploratory | rsf | miRNA only | -0.230337 | [-0.299117, -0.148872] | 0.001998 | 0.003996 |
| exploratory | rsf | Clinical + RNA-seq | -0.149211 | [-0.219107, -0.071048] | 0.001998 | 0.003996 |
| exploratory | rsf | Clinical + RNA-seq + Methylation | -0.142044 | [-0.211959, -0.064209] | 0.001998 | 0.003996 |
| exploratory | rsf | Clinical + RNA-seq + miRNA | -0.144074 | [-0.213223, -0.072581] | 0.001998 | 0.003996 |
| exploratory | rsf | Clinical + Methylation | -0.040110 | [-0.092980, 0.012848] | 0.141858 | 0.148026 |
| exploratory | rsf | Clinical + miRNA | -0.018286 | [-0.048582, 0.014257] | 0.285714 | 0.285714 |
| exploratory | deepsurv | RNA-seq only | -0.172788 | [-0.244014, -0.096515] | 0.001998 | 0.003996 |
| exploratory | deepsurv | Methylation only | -0.288597 | [-0.413680, -0.165460] | 0.001998 | 0.003996 |
| exploratory | deepsurv | miRNA only | -0.245532 | [-0.323178, -0.166095] | 0.001998 | 0.003996 |
| exploratory | deepsurv | Clinical + RNA-seq | -0.130494 | [-0.221154, -0.038912] | 0.011988 | 0.016924 |
| exploratory | deepsurv | Clinical + RNA-seq + Methylation | -0.124250 | [-0.199449, -0.041579] | 0.007992 | 0.013701 |
| exploratory | deepsurv | Clinical + RNA-seq + miRNA | -0.169773 | [-0.262507, -0.080458] | 0.003996 | 0.007377 |
| exploratory | deepsurv | Clinical + Methylation | -0.046231 | [-0.096519, 0.006666] | 0.079920 | 0.091337 |
| exploratory | deepsurv | Clinical + miRNA | -0.159500 | [-0.254934, -0.063264] | 0.001998 | 0.003996 |
