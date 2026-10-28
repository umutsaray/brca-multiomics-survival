# Phase 16 Task 3 — Alternative fusion strategies

Decision rules fixed before results: `{"exploratory": true, "cohort": "corrected primary-tumor-only", "outer_folds": "identical to Table 7", "outer_seed": 42, "selection": "outer-training-only univariate Cox, BH-FDR < 0.05, cap 200/block", "ridge_grid": [0.01, 0.05, 0.1, 0.5, 1.0], "inner_cv": "5-fold participant-grouped/event-stratified; seed=4200+outer_fold", "late_equal": "mean of standardized clinical/RNA/miRNA/methylation Cox risk scores", "late_stack": "ridge CoxPH meta-model on four standardized inner-OOF risks", "pca_k_values": [5, 10, 20], "pca_fit_scope": "outer-training only, separately per omics block", "skipped_runtime": ["RSF-based block late fusion", "PCA DeepSurv"]}`

## C-index

| Configuration | Model | Mean | Sample SD | Fold values |
|---|---|---:|---:|---|
| Clinical-only | coxph | 0.7598208476599109 | 0.07678272079706581 | `[0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375]` |
| Clinical-only | rsf | 0.7572530299391804 | 0.06621912967730519 | `[0.804127198917456, 0.6532967032967033, 0.8178267045454546, 0.7347645429362881, 0.77625]` |
| Full | coxph | 0.6524726851411089 | 0.06777572101415535 | `[0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667]` |
| Full | rsf | 0.6220342521048454 | 0.04251887027710327 | `[0.6427604871447903, 0.6538461538461539, 0.5987215909090909, 0.6565096952908587, 0.5583333333333333]` |
| Late fusion — equal Cox risks | coxph | 0.7092179551332813 | 0.04477419754349126 | `[0.7753721244925575, 0.6582417582417582, 0.7294034090909091, 0.6939058171745153, 0.6891666666666667]` |
| Late fusion — stacked Cox risks | coxph | 0.7028794600559853 | 0.042504672827179145 | `[0.7733423545331529, 0.6648351648351648, 0.7017045454545454, 0.6745152354570637, 0.7]` |
| PCA early fusion k=5 | coxph | 0.7128704242431907 | 0.039345415814810554 | `[0.7571041948579161, 0.7120879120879121, 0.7471590909090909, 0.6696675900277008, 0.6783333333333333]` |
| PCA early fusion k=5 | rsf | 0.6065568593798887 | 0.06307496908250446 | `[0.5656292286874154, 0.6362637362637362, 0.6484375, 0.5166204986149584, 0.6658333333333334]` |
| PCA early fusion k=10 | coxph | 0.6916574720625898 | 0.04440678653353544 | `[0.7368064952638701, 0.7087912087912088, 0.7173295454545454, 0.6703601108033241, 0.625]` |
| PCA early fusion k=10 | rsf | 0.5728098878644573 | 0.058256339755092416 | `[0.5378890392422192, 0.6307692307692307, 0.5852272727272727, 0.4909972299168975, 0.6191666666666666]` |
| PCA early fusion k=20 | coxph | 0.6900609062620214 | 0.04052079969163674 | `[0.7584573748308525, 0.6703296703296703, 0.6953125, 0.6662049861495845, 0.66]` |
| PCA early fusion k=20 | rsf | 0.587749799576811 | 0.05485861285786822 | `[0.5656292286874154, 0.6450549450549451, 0.5653409090909091, 0.5193905817174516, 0.6433333333333333]` |

## Fold-stratified paired bootstrap

| Comparison | Model | Estimate | 95% CI | Raw p | BH q |
|---|---|---:|---|---:|---:|
| Clinical-only minus Late fusion — equal Cox risks | coxph | 0.05060289252662957 | `[-0.010203283596692562, 0.1113524381774762]` | 0.1024 | 0.1664 |
| Late fusion — equal Cox risks minus Full | coxph | 0.0567452699921724 | `[-0.013640213165020985, 0.12730109055490751]` | 0.1144 | 0.1664 |
| Clinical-only minus Late fusion — stacked Cox risks | coxph | 0.05694138760392557 | `[-0.0022499411670819174, 0.11752019991766974]` | 0.0614 | 0.15542857142857144 |
| Late fusion — stacked Cox risks minus Full | coxph | 0.0504067749148764 | `[-0.022102444282566056, 0.11947509901406429]` | 0.1812 | 0.2416 |
| Clinical-only minus PCA early fusion k=5 | coxph | 0.04695042341672027 | `[-0.010757360110012357, 0.10797621233555799]` | 0.112 | 0.1664 |
| Clinical-only minus PCA early fusion k=5 | rsf | 0.15069617055929171 | `[0.07506128901948882, 0.2253714883779904]` | 0.0004 | 0.0021333333333333334 |
| PCA early fusion k=5 minus Full | coxph | 0.0603977391020817 | `[-0.013293814499329844, 0.12982222475631158]` | 0.1062 | 0.1664 |
| PCA early fusion k=5 minus Full | rsf | -0.01547739272495674 | `[-0.06514087717648341, 0.03127002155343887]` | 0.5042 | 0.5042 |
| Clinical-only minus PCA early fusion k=10 | coxph | 0.06816337559732125 | `[0.001954305087243027, 0.13448920754492227]` | 0.0432 | 0.13824 |
| Clinical-only minus PCA early fusion k=10 | rsf | 0.18444314207472304 | `[0.10809056690755176, 0.25738156348476887]` | 0.0002 | 0.0016 |
| PCA early fusion k=10 minus Full | coxph | 0.03918478692148073 | `[-0.03831065523296678, 0.11416530613149116]` | 0.326 | 0.3725714285714286 |
| PCA early fusion k=10 minus Full | rsf | -0.04922436424038806 | `[-0.10242543191477203, 0.003665158267546746]` | 0.068 | 0.15542857142857144 |
| Clinical-only minus PCA early fusion k=20 | coxph | 0.06975994139788948 | `[0.005372648913204143, 0.135658413652871]` | 0.0338 | 0.1352 |
| Clinical-only minus PCA early fusion k=20 | rsf | 0.16950323036236953 | `[0.09466553369251449, 0.24169240888848462]` | 0.0002 | 0.0016 |
| PCA early fusion k=20 minus Full | coxph | 0.03758822112091249 | `[-0.041987077765518785, 0.1147081381347791]` | 0.3638 | 0.38805333333333336 |
| PCA early fusion k=20 minus Full | rsf | -0.03428445252803454 | `[-0.08571127511901235, 0.017576826943935697]` | 0.199 | 0.24492307692307694 |

## Reproducibility

Command: `C:\moenv\Scripts\python.exe src\evaluation\phase16_task3_fusion_strategies.py --input-dir .. --repo-dir . --output-dir results\phase16 --seed 42 --bootstrap-iterations 9999 --meth-prefilter 10000 --knn-neighbors 5 --chunk 5000 --top-n 200 --rsf-trees 500 --max-steps 100`  
Runtime seconds: `1110.4178037999664`.

## CoxPH design-matrix diagnostics

| Fold | Configuration/block | Rows | Columns | Rank | Rank deficient | Convergence warnings |
|---:|---|---:|---:|---:|---|---|
| 1 | Late fusion — equal Cox risks: clinical | 609 | 15 | 12 | True | `[]` |
| 1 | Late fusion — equal Cox risks: rna | 609 | 200 | 191 | True | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 1 | Late fusion — equal Cox risks: mirna | 609 | 199 | 199 | False | `[]` |
| 1 | Late fusion — equal Cox risks: methylation | 609 | 200 | 200 | False | `[]` |
| 1 | Late fusion — stacked Cox risks | 609 | 4 | 4 | False | `[]` |
| 1 | PCA early fusion k=5 | 609 | 30 | 27 | True | `[]` |
| 1 | PCA early fusion k=10 | 609 | 45 | 42 | True | `[]` |
| 1 | PCA early fusion k=20 | 609 | 75 | 72 | True | `[]` |
| 2 | Late fusion — equal Cox risks: clinical | 609 | 14 | 11 | True | `[]` |
| 2 | Late fusion — equal Cox risks: rna | 609 | 200 | 199 | True | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 2 | Late fusion — equal Cox risks: mirna | 609 | 199 | 199 | False | `[]` |
| 2 | Late fusion — equal Cox risks: methylation | 609 | 200 | 200 | False | `[]` |
| 2 | Late fusion — stacked Cox risks | 609 | 4 | 4 | False | `[]` |
| 2 | PCA early fusion k=5 | 609 | 29 | 26 | True | `[]` |
| 2 | PCA early fusion k=10 | 609 | 44 | 41 | True | `[]` |
| 2 | PCA early fusion k=20 | 609 | 74 | 71 | True | `[]` |
| 3 | Late fusion — equal Cox risks: clinical | 610 | 12 | 10 | True | `[]` |
| 3 | Late fusion — equal Cox risks: rna | 610 | 200 | 200 | False | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 3 | Late fusion — equal Cox risks: mirna | 610 | 197 | 197 | False | `[]` |
| 3 | Late fusion — equal Cox risks: methylation | 610 | 200 | 200 | False | `[]` |
| 3 | Late fusion — stacked Cox risks | 610 | 4 | 4 | False | `[]` |
| 3 | PCA early fusion k=5 | 610 | 27 | 25 | True | `[]` |
| 3 | PCA early fusion k=10 | 610 | 42 | 40 | True | `[]` |
| 3 | PCA early fusion k=20 | 610 | 72 | 70 | True | `[]` |
| 4 | Late fusion — equal Cox risks: clinical | 610 | 15 | 12 | True | `[]` |
| 4 | Late fusion — equal Cox risks: rna | 610 | 198 | 196 | True | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 4 | Late fusion — equal Cox risks: mirna | 610 | 198 | 198 | False | `[]` |
| 4 | Late fusion — equal Cox risks: methylation | 610 | 200 | 200 | False | `[]` |
| 4 | Late fusion — stacked Cox risks | 610 | 4 | 4 | False | `[]` |
| 4 | PCA early fusion k=5 | 610 | 30 | 27 | True | `[]` |
| 4 | PCA early fusion k=10 | 610 | 45 | 42 | True | `[]` |
| 4 | PCA early fusion k=20 | 610 | 75 | 72 | True | `[]` |
| 5 | Late fusion — equal Cox risks: clinical | 610 | 15 | 12 | True | `[]` |
| 5 | Late fusion — equal Cox risks: rna | 610 | 200 | 197 | True | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 5 | Late fusion — equal Cox risks: mirna | 610 | 199 | 199 | False | `[]` |
| 5 | Late fusion — equal Cox risks: methylation | 610 | 200 | 200 | False | `[]` |
| 5 | Late fusion — stacked Cox risks | 610 | 4 | 4 | False | `[]` |
| 5 | PCA early fusion k=5 | 610 | 30 | 27 | True | `[]` |
| 5 | PCA early fusion k=10 | 610 | 45 | 42 | True | `[]` |
| 5 | PCA early fusion k=20 | 610 | 75 | 72 | True | `[]` |
