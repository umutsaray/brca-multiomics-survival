# Phase 15 Task 2A — Missing-indicator receptor sensitivity

All decision rules were fixed before model results were inspected. PAM50 was not used as a model input. This is an exploratory sensitivity analysis.

## C-index

| Configuration | Model | Mean | Sample SD | Fold values |
|---|---|---:|---:|---|
| Clinical-only | coxph | 0.7598208476599109 | 0.07678272079706581 | `[0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375]` |
| Clinical-only | rsf | 0.7572530299391804 | 0.06621912967730519 | `[0.804127198917456, 0.6532967032967033, 0.8178267045454546, 0.7347645429362881, 0.77625]` |
| Clinical-only | deepsurv | 0.7443252097099942 | 0.09876928688150217 | `[0.8298376184032477, 0.589010989010989, 0.78515625, 0.7063711911357341, 0.81125]` |
| Clinical + ER/PR/HER2 | coxph | 0.7722063686430005 | 0.068033713365131 | `[0.8328822733423545, 0.6571428571428571, 0.7702414772727273, 0.7995152354570637, 0.80125]` |
| Clinical + ER/PR/HER2 | rsf | 0.8083281712863248 | 0.07627155437637671 | `[0.8585926928281461, 0.6802197802197802, 0.7979403409090909, 0.8653047091412742, 0.8395833333333333]` |
| Clinical + ER/PR/HER2 | deepsurv | 0.690207330919335 | 0.1594280009491826 | `[0.4614343707713126, 0.5923076923076923, 0.7531960227272727, 0.7995152354570637, 0.8445833333333334]` |
| Full | coxph | 0.6524726851411089 | 0.06777572101415535 | `[0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667]` |
| Full | rsf | 0.6220342521048454 | 0.04251887027710327 | `[0.6427604871447903, 0.6538461538461539, 0.5987215909090909, 0.6565096952908587, 0.5583333333333333]` |
| Full | deepsurv | 0.5262759610444295 | 0.03726840910597594 | `[0.5690121786197564, 0.4857142857142857, 0.5426136363636364, 0.5457063711911357, 0.48833333333333334]` |
| Full + ER/PR/HER2 | coxph | 0.6190729620450831 | 0.05395066362393452 | `[0.6657645466847091, 0.554945054945055, 0.6193181818181818, 0.6786703601108033, 0.5766666666666667]` |
| Full + ER/PR/HER2 | rsf | 0.628005753861151 | 0.036355573484750284 | `[0.6515561569688768, 0.6538461538461539, 0.5965909090909091, 0.657202216066482, 0.5808333333333333]` |
| Full + ER/PR/HER2 | deepsurv | 0.5488557112769886 | 0.03556308285291008 | `[0.5967523680649527, 0.5054945054945055, 0.5234375, 0.5560941828254847, 0.5625]` |

## Fold-stratified paired bootstrap

| Comparison | Model | Mean difference | 95% CI | Raw p | BH-FDR q |
|---|---|---:|---|---:|---:|
| Clinical + ER/PR/HER2 minus Clinical-only | coxph | 0.012385520983089603 | `[-0.01349352120942823, 0.04017893523038925]` | 0.3624 | 0.392 |
| Clinical + ER/PR/HER2 minus Clinical-only | rsf | 0.05107514134714457 | `[0.005981696671458864, 0.10034759163799767]` | 0.024 | 0.04114285714285714 |
| Clinical + ER/PR/HER2 minus Clinical-only | deepsurv | -0.05411787879065923 | `[-0.10929741431876268, 0.002152129825225279]` | 0.0606 | 0.09090000000000001 |
| Clinical-only minus Full | coxph | 0.10734816251880197 | `[0.03700289487664033, 0.17268175928129118]` | 0.0034 | 0.0068 |
| Clinical-only minus Full | rsf | 0.13521877783433497 | `[0.05837150819325333, 0.20921014818192593]` | 0.0008 | 0.00192 |
| Clinical-only minus Full | deepsurv | 0.21804924866556466 | `[0.11942560134539079, 0.3050861324581443]` | 0.0002 | 0.0008 |
| Clinical + ER/PR/HER2 minus Full | coxph | 0.11973368350189158 | `[0.05464276440556956, 0.18160474010669803]` | 0.0002 | 0.0008 |
| Clinical + ER/PR/HER2 minus Full | rsf | 0.18629391918147956 | `[0.11687171578464438, 0.25228298307241004]` | 0.0002 | 0.0008 |
| Clinical + ER/PR/HER2 minus Full | deepsurv | 0.16393136987490542 | `[0.06733884616167264, 0.2503614375310016]` | 0.0006 | 0.0018 |
| Full + ER/PR/HER2 minus Full | coxph | -0.03339972309602581 | `[-0.07175286424850096, 0.010246461851344514]` | 0.1288 | 0.17173333333333332 |
| Full + ER/PR/HER2 minus Full | rsf | 0.0059715017563055865 | `[-0.008458316557561716, 0.02062921199537629]` | 0.392 | 0.392 |
| Full + ER/PR/HER2 minus Full | deepsurv | 0.022579750232559072 | `[-0.025699732098941386, 0.06957833329470663]` | 0.3712 | 0.392 |

## Design-matrix and convergence diagnostics

| Fold | Configuration | Model | Rows | Columns | Rank | Rank deficient | Convergence warnings |
|---:|---|---|---:|---:|---:|---|---|
| 1 | Clinical-only | coxph | 609 | 15 | 12 | true | `[]` |
| 1 | Clinical-only | rsf | 609 | 16 | 13 | true | `[]` |
| 1 | Clinical-only | deepsurv | 522 | 16 | 13 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | coxph | 609 | 21 | 18 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | rsf | 609 | 22 | 19 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | deepsurv | 522 | 22 | 19 | true | `[]` |
| 1 | Full | coxph | 609 | 217 | 205 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 1 | Full | rsf | 609 | 218 | 191 | true | `[]` |
| 1 | Full | deepsurv | 522 | 229 | 202 | true | `[]` |
| 1 | Full + ER/PR/HER2 | coxph | 609 | 223 | 211 | true | `[]` |
| 1 | Full + ER/PR/HER2 | rsf | 609 | 224 | 197 | true | `[]` |
| 1 | Full + ER/PR/HER2 | deepsurv | 522 | 235 | 208 | true | `[]` |
| 2 | Clinical-only | coxph | 609 | 14 | 11 | true | `[]` |
| 2 | Clinical-only | rsf | 609 | 15 | 12 | true | `[]` |
| 2 | Clinical-only | deepsurv | 522 | 15 | 12 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | coxph | 609 | 20 | 17 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | rsf | 609 | 21 | 18 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | deepsurv | 522 | 21 | 18 | true | `[]` |
| 2 | Full | coxph | 609 | 220 | 216 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 2 | Full | rsf | 609 | 221 | 207 | true | `[]` |
| 2 | Full | deepsurv | 522 | 230 | 199 | true | `[]` |
| 2 | Full + ER/PR/HER2 | coxph | 609 | 226 | 222 | true | `[]` |
| 2 | Full + ER/PR/HER2 | rsf | 609 | 227 | 213 | true | `[]` |
| 2 | Full + ER/PR/HER2 | deepsurv | 522 | 236 | 205 | true | `[]` |
| 3 | Clinical-only | coxph | 610 | 12 | 10 | true | `[]` |
| 3 | Clinical-only | rsf | 610 | 16 | 13 | true | `[]` |
| 3 | Clinical-only | deepsurv | 523 | 16 | 13 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | coxph | 610 | 18 | 16 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | rsf | 610 | 22 | 19 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | deepsurv | 523 | 22 | 19 | true | `[]` |
| 3 | Full | coxph | 610 | 140 | 138 | true | `[]` |
| 3 | Full | rsf | 610 | 145 | 138 | true | `[]` |
| 3 | Full | deepsurv | 523 | 159 | 153 | true | `[]` |
| 3 | Full + ER/PR/HER2 | coxph | 610 | 146 | 144 | true | `[]` |
| 3 | Full + ER/PR/HER2 | rsf | 610 | 151 | 144 | true | `[]` |
| 3 | Full + ER/PR/HER2 | deepsurv | 523 | 165 | 159 | true | `[]` |
| 4 | Clinical-only | coxph | 610 | 15 | 12 | true | `[]` |
| 4 | Clinical-only | rsf | 610 | 16 | 13 | true | `[]` |
| 4 | Clinical-only | deepsurv | 522 | 16 | 13 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | coxph | 610 | 21 | 18 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | rsf | 610 | 22 | 19 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | deepsurv | 522 | 22 | 19 | true | `[]` |
| 4 | Full | coxph | 610 | 227 | 221 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 4 | Full | rsf | 610 | 230 | 199 | true | `[]` |
| 4 | Full | deepsurv | 522 | 233 | 200 | true | `[]` |
| 4 | Full + ER/PR/HER2 | coxph | 610 | 233 | 227 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 4 | Full + ER/PR/HER2 | rsf | 610 | 236 | 205 | true | `[]` |
| 4 | Full + ER/PR/HER2 | deepsurv | 522 | 239 | 206 | true | `[]` |
| 5 | Clinical-only | coxph | 610 | 15 | 12 | true | `[]` |
| 5 | Clinical-only | rsf | 610 | 16 | 13 | true | `[]` |
| 5 | Clinical-only | deepsurv | 522 | 16 | 13 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | coxph | 610 | 21 | 18 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | rsf | 610 | 22 | 19 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | deepsurv | 522 | 22 | 19 | true | `[]` |
| 5 | Full | coxph | 610 | 231 | 225 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 5 | Full | rsf | 610 | 232 | 212 | true | `[]` |
| 5 | Full | deepsurv | 522 | 135 | 127 | true | `[]` |
| 5 | Full + ER/PR/HER2 | coxph | 610 | 237 | 231 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 5 | Full + ER/PR/HER2 | rsf | 610 | 238 | 218 | true | `[]` |
| 5 | Full + ER/PR/HER2 | deepsurv | 522 | 141 | 133 | true | `[]` |

## Table 7 reproduction check

`{"Clinical-only": {"coxph": {"observed": [0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375], "table7": [0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}, "rsf": {"observed": [0.804127198917456, 0.6532967032967033, 0.8178267045454546, 0.7347645429362881, 0.77625], "table7": [0.804127198917456, 0.6532967032967033, 0.8178267045454546, 0.7347645429362881, 0.77625], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}, "deepsurv": {"observed": [0.8298376184032477, 0.589010989010989, 0.78515625, 0.7063711911357341, 0.81125], "table7": [0.8298376184032477, 0.589010989010989, 0.78515625, 0.7063711911357341, 0.81125], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}}, "Full": {"coxph": {"observed": [0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667], "table7": [0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}, "rsf": {"observed": [0.6427604871447903, 0.6538461538461539, 0.5987215909090909, 0.6565096952908587, 0.5583333333333333], "table7": [0.6427604871447903, 0.6538461538461539, 0.5987215909090909, 0.6565096952908587, 0.5583333333333333], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}, "deepsurv": {"observed": [0.5690121786197564, 0.4857142857142857, 0.5426136363636364, 0.5457063711911357, 0.48833333333333334], "table7": [0.5690121786197564, 0.4857142857142857, 0.5426136363636364, 0.5457063711911357, 0.48833333333333334], "exact_match": true, "differences": [0.0, 0.0, 0.0, 0.0, 0.0]}}}`

## Reproducibility

Command: `C:\moenv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase15_task2_receptor_sensitivity.py --input-dir . --repo-dir .\brca-multiomics-survival --output-dir .\brca-multiomics-survival\results\phase15 --outer-seed 42 --bootstrap-iterations 9999 --bootstrap-seed 42 --meth-prefilter 10000 --knn-neighbors 5 --methylation-chunksize 5000 --top-n 200 --cox-penalizer 0.1 --cox-max-steps 50 --rsf-trees 500 --dl-epochs 300 --dl-patience 30 --dl-learning-rate 0.001 --dl-weight-decay 0.0001 --dl-dropout 0.3`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3", "numpy": "2.5.2", "scikit_learn": "1.9.0", "scikit_survival": "0.28.0", "lifelines": "0.30.3", "statsmodels": "0.15.0", "torch": "2.13.0+cpu"}`  
Seeds: `{"outer": 42, "deepsurv_internal_schema": "4200 + outer_fold", "model_schema": "42 + outer_fold", "bootstrap": 42}`  
Runtime seconds: `638.6391574000008`.
