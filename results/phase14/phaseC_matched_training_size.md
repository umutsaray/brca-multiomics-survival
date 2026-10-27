# Phase 14 Check C — Matched training-set size

Status: **completed**.

The exact Phase 7 internal fit subsets and cached preprocessing/feature selection were reused. Phase 7 holds out one of seven grouped folds (~14.3%), so the retained fit subset is described as matched-85% but is not a newly generated exact 15% holdout.

## C-index results

| Configuration | Model | Full outer-training mean | Full SD | Matched fit-subset mean | Matched SD | Matched fold C-index values |
|---|---|---:|---:|---:|---:|---|
| Clinical-only | coxph | 0.7598208476599109 | 0.07678272079706581 | 0.772009163957373 | 0.0815888213992139 | `[0.8406630581867388, 0.6384615384615384, 0.8263494318181818, 0.7541551246537396, 0.8004166666666667]` |
| Clinical-only | rsf | 0.7572530299391804 | 0.06621912967730519 | 0.7534373254009958 | 0.07156129509386909 | `[0.8146143437077131, 0.6412087912087913, 0.8121448863636364, 0.7313019390581718, 0.7679166666666667]` |
| Full | coxph | 0.6524726851411089 | 0.06777572101415535 | 0.6474605487882941 | 0.04798691956923946 | `[0.6684709066305818, 0.5769230769230769, 0.6640625, 0.6253462603878116, 0.7025]` |
| Full | rsf | 0.6220342521048454 | 0.04251887027710327 | 0.6260923017409405 | 0.08300266995857804 | `[0.699594046008119, 0.5494505494505495, 0.5262784090909091, 0.7001385041551247, 0.655]` |

## Fold-stratified paired bootstrap: matched minus full

| Configuration | Model | Mean difference | 95% CI | Two-sided p | Iterations | Seed |
|---|---|---:|---|---:|---:|---:|
| Clinical-only | coxph | 0.01218831629746211 | `[-0.0019277696323236778, 0.027766595870547948]` | 0.0914 | 9999 | 42 |
| Clinical-only | rsf | -0.003815704538184561 | `[-0.020785119436129953, 0.012682951771289933]` | 0.6906 | 9999 | 42 |
| Full | coxph | -0.0050121363528149 | `[-0.08265230568839338, 0.07277360580288222]` | 0.906 | 9999 | 42 |
| Full | rsf | 0.0040580496360950406 | `[-0.047266189309712575, 0.052696575075210215]` | 0.8696 | 9999 | 42 |

## Ordering checks

Clinical-only versus Full:

- `coxph`: changed = `false`; full-training Clinical−Full = `0.10734816251880197`; matched Clinical−Full = `0.12454861516907889`.
- `rsf`: changed = `false`; full-training Clinical−Full = `0.135218777834335`; matched Clinical−Full = `0.1273450236600553`.

DeepSurv versus classical model:

- `Clinical-only`, DeepSurv vs `coxph`: changed = `false`; before DeepSurv−classical = `-0.015495637949916707`; after = `-0.027683954247378795`.
- `Clinical-only`, DeepSurv vs `rsf`: changed = `false`; before DeepSurv−classical = `-0.01292782022918626`; after = `-0.009112115691001677`.
- `Full`, DeepSurv vs `coxph`: changed = `false`; before DeepSurv−classical = `-0.12619672409667937`; after = `-0.12118458774386454`.
- `Full`, DeepSurv vs `rsf`: changed = `false`; before DeepSurv−classical = `-0.09575829106041589`; after = `-0.09981634069651102`.

## Fold details

- Fold `1`: outer train `609`, fit subset `522`, internal validation `87`, outer test `153`, internal seed `4201`, model seed `43`.
- Fold `2`: outer train `609`, fit subset `522`, internal validation `87`, outer test `153`, internal seed `4202`, model seed `44`.
- Fold `3`: outer train `610`, fit subset `523`, internal validation `87`, outer test `152`, internal seed `4203`, model seed `45`.
- Fold `4`: outer train `610`, fit subset `522`, internal validation `88`, outer test `152`, internal seed `4204`, model seed `46`.
- Fold `5`: outer train `610`, fit subset `522`, internal validation `88`, outer test `152`, internal seed `4205`, model seed `47`.

## Reproducibility

Command: `C:\moenv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase14_c_matched_training_size.py --input-dir . --repo-dir .\brca-multiomics-survival --output-dir .\brca-multiomics-survival\results\phase14 --outer-seed 42 --bootstrap-iterations 9999 --bootstrap-seed 42 --cox-penalizer 0.1 --cox-max-steps 50 --rsf-trees 500`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3", "numpy": "2.5.2", "scikit_learn": "1.9.0", "scikit_survival": "0.28.0", "lifelines": "0.30.3"}`  
Runtime seconds: `161.79984889994375`.
