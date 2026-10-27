# Phase 14 quick diagnostic summary

## Check A — ER / PR / HER2 availability

Status: **completed**.
ER: available `458`, missing `304`, events with value `64` / `99`; categories `{"Negative": 103, "Positive": 355}`.
PR: available `457`, missing `305`, events with value `64` / `99`; categories `{"Indeterminate": 2, "Negative": 140, "Positive": 315}`.
HER2: available `455`, missing `307`, events with value `59` / `99`; categories `{"Equivocal": 3, "Negative": 395, "Positive": 57}`.
PAM50: available `523`, missing `239`, events with value `73` / `99`; categories `{"Basal": 78, "Her2": 29, "LumA": 275, "LumB": 124, "Normal": 17}`.

## Check B — RNA-seq scale

Status: **stopped_no_modelling_needed**.
The matrix is already log2(TPM+1), as defined by the UCSC Xena GDC STAR-TPM ETL; no additional log transform was applied.
Observed minimum `0.0`, maximum `18.78925718376024`, median `0.0716253939737696`; percentiles `{"0": 0.0, "1": 0.0, "5": 0.0, "10": 0.0, "25": 0.0, "50": 0.0716253939737696, "75": 1.5191392826429724, "90": 4.512770481607714, "95": 5.6619291275694135, "99": 7.686658094493024, "100": 18.78925718376024}`.

## Check C — Matched training size

Status: **completed**.
| Configuration | Model | Full mean ± SD | Matched mean ± SD | Matched−full mean difference | 95% CI | p |
|---|---|---|---|---:|---|---:|
| Clinical-only | coxph | `0.7598208476599109 ± 0.07678272079706581` | `0.772009163957373 ± 0.0815888213992139` | 0.01218831629746211 | `[-0.0019277696323236778, 0.027766595870547948]` | 0.0914 |
| Clinical-only | rsf | `0.7572530299391804 ± 0.06621912967730519` | `0.7534373254009958 ± 0.07156129509386909` | -0.003815704538184561 | `[-0.020785119436129953, 0.012682951771289933]` | 0.6906 |
| Full | coxph | `0.6524726851411089 ± 0.06777572101415535` | `0.6474605487882941 ± 0.04798691956923946` | -0.0050121363528149 | `[-0.08265230568839338, 0.07277360580288222]` | 0.906 |
| Full | rsf | `0.6220342521048454 ± 0.04251887027710327` | `0.6260923017409405 ± 0.08300266995857804` | 0.0040580496360950406 | `[-0.047266189309712575, 0.052696575075210215]` | 0.8696 |

Ordering changes:

- Clinical-only vs Full, `coxph`: `false` (full-training delta `0.10734816251880197`; matched delta `0.12454861516907889`).
- Clinical-only vs Full, `rsf`: `false` (full-training delta `0.135218777834335`; matched delta `0.1273450236600553`).
- `Clinical-only`, DeepSurv vs `coxph`: `false` (before `-0.015495637949916707`; after `-0.027683954247378795`).
- `Clinical-only`, DeepSurv vs `rsf`: `false` (before `-0.01292782022918626`; after `-0.009112115691001677`).
- `Full`, DeepSurv vs `coxph`: `false` (before `-0.12619672409667937`; after `-0.12118458774386454`).
- `Full`, DeepSurv vs `rsf`: `false` (before `-0.09575829106041589`; after `-0.09981634069651102`).

## Reproducibility records

- Check A command: `C:\Users\UMUT\Desktop\Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer\.venv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase14_a_clinical_availability.py --input-dir . --output-dir .\brca-multiomics-survival\results\phase14`; runtime seconds `2.6271913999225944`; seeds `null`.
- Check B command: `C:\Users\UMUT\Desktop\Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer\.venv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase14_b_rnaseq_scale.py --input-dir . --output-dir .\brca-multiomics-survival\results\phase14`; runtime seconds `13.943086200160906`; seeds `null`.
- Check C command: `C:\moenv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase14_c_matched_training_size.py --input-dir . --repo-dir .\brca-multiomics-survival --output-dir .\brca-multiomics-survival\results\phase14 --outer-seed 42 --bootstrap-iterations 9999 --bootstrap-seed 42 --cox-penalizer 0.1 --cox-max-steps 50 --rsf-trees 500`; runtime seconds `161.79984889994375`; outer seed `42`; bootstrap seed `42`.
