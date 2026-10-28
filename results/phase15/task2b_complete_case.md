# Phase 15 Task 2B — Complete-case receptor sensitivity

All decision rules were fixed before model results were inspected. PAM50 was not used as a model input. This is an exploratory sensitivity analysis.

## C-index

| Configuration | Model | Mean | Sample SD | Fold values |
|---|---|---:|---:|---|
| Clinical-only | coxph | 0.7368866149723537 | 0.05158508321234812 | `[0.691304347826087, 0.696078431372549, 0.813953488372093, 0.7201210287443268, 0.7629757785467128]` |
| Clinical-only | rsf | 0.7077419855122589 | 0.09988233007840161 | `[0.6130434782608696, 0.5857843137254902, 0.7921511627906976, 0.7639939485627837, 0.7837370242214533]` |
| Clinical-only | deepsurv | 0.6470473735705957 | 0.06258835968847383 | `[0.591304347826087, 0.5980392156862745, 0.6555232558139535, 0.642965204236006, 0.7474048442906575]` |
| Clinical + ER/PR/HER2 | coxph | 0.7396397606257427 | 0.06188574902617041 | `[0.6565217391304348, 0.7640931372549019, 0.8241279069767442, 0.7397881996974282, 0.7136678200692042]` |
| Clinical + ER/PR/HER2 | rsf | 0.7562771683057392 | 0.07498210259512648 | `[0.6652173913043479, 0.7628676470588235, 0.8619186046511628, 0.7829046898638427, 0.7084775086505191]` |
| Clinical + ER/PR/HER2 | deepsurv | 0.7081131223014553 | 0.046435266066052534 | `[0.6282608695652174, 0.7469362745098039, 0.7194767441860465, 0.7322239031770046, 0.7136678200692042]` |
| Full | coxph | 0.6055094029868497 | 0.11895712990775721 | `[0.758695652173913, 0.4375, 0.6671511627906976, 0.5794251134644478, 0.5847750865051903]` |
| Full | rsf | 0.6080089585487123 | 0.073276100985532 | `[0.6869565217391305, 0.6715686274509803, 0.6148255813953488, 0.5355521936459909, 0.5311418685121108]` |
| Full | deepsurv | 0.5511203940914301 | 0.10400081808829793 | `[0.5891304347826087, 0.5772058823529411, 0.6773255813953488, 0.39636913767019666, 0.5155709342560554]` |
| Full + ER/PR/HER2 | coxph | 0.6346160389320257 | 0.06998528090883313 | `[0.7434782608695653, 0.5906862745098039, 0.6671511627906976, 0.5869894099848714, 0.5847750865051903]` |
| Full + ER/PR/HER2 | rsf | 0.6009000648110427 | 0.06934291011640975 | `[0.6652173913043479, 0.6580882352941176, 0.627906976744186, 0.5325264750378215, 0.5207612456747405]` |
| Full + ER/PR/HER2 | deepsurv | 0.5121643992179209 | 0.08987515988596162 | `[0.5934782608695652, 0.5477941176470589, 0.5334302325581395, 0.5279878971255674, 0.35813148788927335]` |

## Fold-stratified paired bootstrap

| Comparison | Model | Mean difference | 95% CI | Raw p | BH-FDR q |
|---|---|---:|---|---:|---:|
| Clinical + ER/PR/HER2 minus Clinical-only | coxph | 0.002753145653388911 | `[-0.03309288917110469, 0.039489017609158315]` | 0.8932 | 0.8932 |
| Clinical + ER/PR/HER2 minus Clinical-only | rsf | 0.04853518279348028 | `[-0.004049391236089664, 0.10226253176596002]` | 0.0728 | 0.1456 |
| Clinical + ER/PR/HER2 minus Clinical-only | deepsurv | 0.06106574873085961 | `[-0.03933864882284523, 0.15095879746449536]` | 0.232 | 0.34800000000000003 |
| Clinical-only minus Full | coxph | 0.13137721198550398 | `[0.030340648901637145, 0.2247032680354341]` | 0.0124 | 0.0496 |
| Clinical-only minus Full | rsf | 0.09973302696354662 | `[-0.006458688769832033, 0.2009666481735296]` | 0.0634 | 0.1456 |
| Clinical-only minus Full | deepsurv | 0.09592697947916555 | `[-0.03964436468159284, 0.21903097893155005]` | 0.159 | 0.2725714285714286 |
| Clinical + ER/PR/HER2 minus Full | coxph | 0.1341303576388929 | `[0.03625885711439692, 0.2227281353098604]` | 0.008 | 0.048 |
| Clinical + ER/PR/HER2 minus Full | rsf | 0.1482682097570269 | `[0.06106784793853994, 0.23005854947295631]` | 0.0016 | 0.019200000000000002 |
| Clinical + ER/PR/HER2 minus Full | deepsurv | 0.15699272821002516 | `[0.020311855888499583, 0.2686434190406548]` | 0.0242 | 0.0726 |
| Full + ER/PR/HER2 minus Full | coxph | 0.02910663594517595 | `[-0.025762781965245788, 0.08341730236370788]` | 0.3148 | 0.37776 |
| Full + ER/PR/HER2 minus Full | rsf | -0.007108893737669564 | `[-0.027809732670908135, 0.013289480749432522]` | 0.499 | 0.5443636363636364 |
| Full + ER/PR/HER2 minus Full | deepsurv | -0.03895599487350928 | `[-0.11093665921259438, 0.02692144370933518]` | 0.2664 | 0.3552 |

## Design-matrix and convergence diagnostics

| Fold | Configuration | Model | Rows | Columns | Rank | Rank deficient | Convergence warnings |
|---:|---|---|---:|---:|---:|---|---|
| 1 | Clinical-only | coxph | 424 | 13 | 11 | true | `[]` |
| 1 | Clinical-only | rsf | 424 | 16 | 13 | true | `[]` |
| 1 | Clinical-only | deepsurv | 363 | 16 | 13 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | coxph | 424 | 16 | 14 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | rsf | 424 | 22 | 16 | true | `[]` |
| 1 | Clinical + ER/PR/HER2 | deepsurv | 363 | 22 | 16 | true | `[]` |
| 1 | Full | coxph | 424 | 181 | 179 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 1 | Full | rsf | 424 | 190 | 171 | true | `[]` |
| 1 | Full | deepsurv | 363 | 211 | 188 | true | `[]` |
| 1 | Full + ER/PR/HER2 | coxph | 424 | 184 | 182 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 1 | Full + ER/PR/HER2 | rsf | 424 | 196 | 174 | true | `[]` |
| 1 | Full + ER/PR/HER2 | deepsurv | 363 | 217 | 191 | true | `[]` |
| 2 | Clinical-only | coxph | 424 | 11 | 10 | true | `[]` |
| 2 | Clinical-only | rsf | 424 | 16 | 13 | true | `[]` |
| 2 | Clinical-only | deepsurv | 363 | 16 | 13 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | coxph | 424 | 14 | 13 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | rsf | 424 | 22 | 16 | true | `[]` |
| 2 | Clinical + ER/PR/HER2 | deepsurv | 363 | 22 | 16 | true | `[]` |
| 2 | Full | coxph | 424 | 140 | 139 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 2 | Full | rsf | 424 | 149 | 134 | true | `[]` |
| 2 | Full | deepsurv | 363 | 148 | 137 | true | `[]` |
| 2 | Full + ER/PR/HER2 | coxph | 424 | 143 | 142 | true | `[]` |
| 2 | Full + ER/PR/HER2 | rsf | 424 | 155 | 137 | true | `[]` |
| 2 | Full + ER/PR/HER2 | deepsurv | 363 | 154 | 140 | true | `[]` |
| 3 | Clinical-only | coxph | 424 | 13 | 11 | true | `[]` |
| 3 | Clinical-only | rsf | 424 | 16 | 13 | true | `[]` |
| 3 | Clinical-only | deepsurv | 363 | 16 | 13 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | coxph | 424 | 16 | 14 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | rsf | 424 | 22 | 16 | true | `[]` |
| 3 | Clinical + ER/PR/HER2 | deepsurv | 363 | 22 | 16 | true | `[]` |
| 3 | Full | coxph | 424 | 161 | 159 | true | `[]` |
| 3 | Full | rsf | 424 | 165 | 156 | true | `[]` |
| 3 | Full | deepsurv | 363 | 206 | 191 | true | `[]` |
| 3 | Full + ER/PR/HER2 | coxph | 424 | 164 | 162 | true | `[]` |
| 3 | Full + ER/PR/HER2 | rsf | 424 | 171 | 159 | true | `[]` |
| 3 | Full + ER/PR/HER2 | deepsurv | 363 | 212 | 194 | true | `[]` |
| 4 | Clinical-only | coxph | 424 | 12 | 11 | true | `[]` |
| 4 | Clinical-only | rsf | 424 | 16 | 13 | true | `[]` |
| 4 | Clinical-only | deepsurv | 363 | 16 | 13 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | coxph | 424 | 15 | 14 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | rsf | 424 | 22 | 16 | true | `[]` |
| 4 | Clinical + ER/PR/HER2 | deepsurv | 363 | 22 | 16 | true | `[]` |
| 4 | Full | coxph | 424 | 182 | 181 | true | `[]` |
| 4 | Full | rsf | 424 | 192 | 176 | true | `[]` |
| 4 | Full | deepsurv | 363 | 229 | 209 | true | `[]` |
| 4 | Full + ER/PR/HER2 | coxph | 424 | 185 | 184 | true | `[]` |
| 4 | Full + ER/PR/HER2 | rsf | 424 | 198 | 178 | true | `[]` |
| 4 | Full + ER/PR/HER2 | deepsurv | 363 | 235 | 212 | true | `[]` |
| 5 | Clinical-only | coxph | 424 | 13 | 10 | true | `[]` |
| 5 | Clinical-only | rsf | 424 | 15 | 12 | true | `[]` |
| 5 | Clinical-only | deepsurv | 363 | 15 | 12 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | coxph | 424 | 16 | 13 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | rsf | 424 | 21 | 15 | true | `[]` |
| 5 | Clinical + ER/PR/HER2 | deepsurv | 363 | 21 | 15 | true | `[]` |
| 5 | Full | coxph | 424 | 198 | 195 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 5 | Full | rsf | 424 | 209 | 188 | true | `[]` |
| 5 | Full | deepsurv | 363 | 209 | 191 | true | `[]` |
| 5 | Full + ER/PR/HER2 | coxph | 424 | 201 | 198 | true | `["Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model"]` |
| 5 | Full + ER/PR/HER2 | rsf | 424 | 215 | 191 | true | `[]` |
| 5 | Full + ER/PR/HER2 | deepsurv | 363 | 215 | 194 | true | `[]` |

## Complete-case fold event counts

Counts: `[14, 15, 15, 14, 14]`.
Folds with fewer than five test events: `[]`.

## Reproducibility

Command: `C:\moenv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase15_task2_receptor_sensitivity.py --input-dir . --repo-dir .\brca-multiomics-survival --output-dir .\brca-multiomics-survival\results\phase15 --outer-seed 42 --bootstrap-iterations 9999 --bootstrap-seed 42 --meth-prefilter 10000 --knn-neighbors 5 --methylation-chunksize 5000 --top-n 200 --cox-penalizer 0.1 --cox-max-steps 50 --rsf-trees 500 --dl-epochs 300 --dl-patience 30 --dl-learning-rate 0.001 --dl-weight-decay 0.0001 --dl-dropout 0.3`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3", "numpy": "2.5.2", "scikit_learn": "1.9.0", "scikit_survival": "0.28.0", "lifelines": "0.30.3", "statsmodels": "0.15.0", "torch": "2.13.0+cpu"}`  
Seeds: `{"outer": 42, "deepsurv_internal_schema": "4200 + outer_fold", "model_schema": "42 + outer_fold", "bootstrap": 42}`  
Runtime seconds: `6310.380588200176`.
