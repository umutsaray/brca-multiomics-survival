# DeepSurv internal-validation leakage correction

Phase 7 splits each outer-training set into participant-grouped, event-stratified internal-train and internal-validation subsets before any learned preprocessing or outcome-guided feature selection. Internal validation is used only for early stopping; each outer test fold is evaluated once.

Scope: exactly 20 retrainings (2 cohorts × Clinical-only/Full × 5 folds). Other ablation configurations were not rerun.

## Old versus new DeepSurv C-index

| Cohort | Configuration | Fold | Old C | New C | New − old |
|---|---|---:|---:|---:|---:|
| mixed | Clinical-only | 1 | 0.6852 | 0.6833 | -0.0019 |
| mixed | Clinical-only | 2 | 0.6967 | 0.7052 | 0.0086 |
| mixed | Clinical-only | 3 | 0.4681 | 0.7362 | 0.2681 |
| mixed | Clinical-only | 4 | 0.7152 | 0.7811 | 0.0659 |
| mixed | Clinical-only | 5 | 0.7883 | 0.7917 | 0.0035 |
| mixed | Full | 1 | 0.5268 | 0.5790 | 0.0521 |
| mixed | Full | 2 | 0.6956 | 0.6994 | 0.0039 |
| mixed | Full | 3 | 0.6022 | 0.4326 | -0.1696 |
| mixed | Full | 4 | 0.6213 | 0.6303 | 0.0090 |
| mixed | Full | 5 | 0.6594 | 0.6191 | -0.0403 |
| primary | Clinical-only | 1 | 0.8248 | 0.8298 | 0.0051 |
| primary | Clinical-only | 2 | 0.5962 | 0.5890 | -0.0071 |
| primary | Clinical-only | 3 | 0.8161 | 0.7852 | -0.0309 |
| primary | Clinical-only | 4 | 0.7140 | 0.7064 | -0.0076 |
| primary | Clinical-only | 5 | 0.7937 | 0.8113 | 0.0175 |
| primary | Full | 1 | 0.5947 | 0.5690 | -0.0257 |
| primary | Full | 2 | 0.5538 | 0.4857 | -0.0681 |
| primary | Full | 3 | 0.4702 | 0.5426 | 0.0724 |
| primary | Full | 4 | 0.5485 | 0.5457 | -0.0028 |
| primary | Full | 5 | 0.4967 | 0.4883 | -0.0083 |

Mean ± sample SD:

| Cohort | Configuration | Old | New |
|---|---|---:|---:|
| mixed | Clinical-only | 0.6707 ± 0.1201 | 0.7395 ± 0.0469 |
| mixed | Full | 0.6211 ± 0.0637 | 0.5921 ± 0.0991 |
| primary | Clinical-only | 0.7489 ± 0.0960 | 0.7443 ± 0.0988 |
| primary | Full | 0.5328 ± 0.0494 | 0.5263 ± 0.0373 |

## SHAP old versus new

Old cross-fold mean pairwise top-20 Jaccard: `0.023627676259255204`.
New cross-fold mean pairwise top-20 Jaccard: `0.028890834153992052`.
Old/new aggregate top-20 overlap: `1/20`; Jaccard: `0.02564102564102564`.

| Rank | Old aggregate feature | New aggregate feature | New mean |SHAP| | New folds selected |
|---:|---|---|---:|---:|
| 1 | rna_ENSG00000233304.7 | meth_cg16596052 | 0.92206056 | 1/5 |
| 2 | meth_cg10973622 | rna_ENSG00000169139.12 | 0.54729591 | 1/5 |
| 3 | meth_cg22001496 | meth_cg06430465 | 0.50038042 | 1/5 |
| 4 | meth_cg04985582 | rna_ENSG00000224177.8 | 0.27405734 | 1/5 |
| 5 | mirna_hsa-mir-7108 | rna_ENSG00000137547.9 | 0.25587504 | 1/5 |
| 6 | mirna_hsa-mir-302b | meth_cg16362594 | 0.22630409 | 1/5 |
| 7 | rna_ENSG00000250765.6 | mirna_hsa-mir-548an | 0.214361 | 1/5 |
| 8 | rna_ENSG00000228361.2 | meth_cg12462916 | 0.20863951 | 1/5 |
| 9 | rna_ENSG00000122025.15 | rna_ENSG00000235209.1 | 0.20200006 | 1/5 |
| 10 | rna_ENSG00000253105.6 | rna_ENSG00000244694.7 | 0.18840211 | 1/5 |
| 11 | meth_cg23866687 | rna_ENSG00000234089.2 | 0.18126385 | 1/5 |
| 12 | rna_ENSG00000072080.11 | rna_ENSG00000227748.2 | 0.17885516 | 2/5 |
| 13 | rna_ENSG00000285715.1 | rna_ENSG00000228361.2 | 0.1770508 | 1/5 |
| 14 | rna_ENSG00000183255.12 | rna_ENSG00000102144.15 | 0.17184399 | 2/5 |
| 15 | mirna_hsa-mir-7107 | rna_ENSG00000147509.14 | 0.16459225 | 1/5 |
| 16 | meth_cg07972135 | rna_ENSG00000235893.5 | 0.16330022 | 1/5 |
| 17 | rna_ENSG00000229145.1 | rna_ENSG00000267484.1 | 0.16113933 | 1/5 |
| 18 | meth_cg04747226 | rna_ENSG00000228974.5 | 0.16057113 | 1/5 |
| 19 | mirna_hsa-mir-6846 | rna_ENSG00000288013.1 | 0.15918982 | 1/5 |
| 20 | rna_ENSG00000076043.11 | rna_ENSG00000231894.2 | 0.15572515 | 1/5 |

## Primary finding after correction

- mixed: Clinical − Full = 0.1474 (95% CI [0.0518, 0.2388], p=0.0034); Clinical-only superiority retained.
- primary: Clinical − Full = 0.2180 (95% CI [0.1194, 0.3051], p=0.0002); Clinical-only superiority retained.

## Seeds and split schema

- Outer StratifiedGroupKFold: seed 42.
- Internal 7-fold StratifiedGroupKFold first held split: seed = 4200 + outer-fold (4201–4205).
- DeepSurv initialization/training: seed = 42 + outer-fold (43–47).
- Fold-stratified paired participant bootstrap: seed 42, 9,999 iterations.

CoxPH and RSF were not retrained. Their Table 9 metrics were recomputed from retained OOF predictions and matched the prior results with maximum absolute difference 0.0.