# Phase 15 Task 1 — ER/PR/HER2 missingness

## Pre-specified decision rules (fixed before model results)

- Source order: `{"ER": ["breast_carcinoma_estrogen_receptor_status", "ER_Status_nature2012"], "PR": ["breast_carcinoma_progesterone_receptor_status", "PR_Status_nature2012"], "HER2": ["lab_proc_her2_neu_immunohistochemistry_receptor_status", "lab_procedure_her2_neu_in_situ_hybrid_outcome_type", "HER2_Final_Status_nature2012"]}`.
- Use the first nonblank source only when every earlier source is blank.
- If available sources contain different normalized nonblank categories, mark `Conflict`; do not resolve it.
- Only `Positive` and `Negative` are model-available. Every other result is `Not available`.
- PAM50 is reported only and is not a model input.

## Candidate columns and counts in the 762-participant cohort

| Column | Non-missing | Blank/NA | Multiple distinct values | Full value counts |
|---|---:|---:|---:|---|
| `ER_Status_nature2012` | 458 | 304 | 0 | `{"Negative": 103, "Positive": 355, "blank/NA": 304}` |
| `HER2_Final_Status_nature2012` | 455 | 307 | 0 | `{"Equivocal": 3, "Negative": 395, "Positive": 57, "blank/NA": 307}` |
| `Integrated_Clusters_with_PAM50__nature2012` | 140 | 622 | 0 | `{"1.0": 7, "2.0": 32, "3.0": 72, "4.0": 29, "blank/NA": 622}` |
| `PAM50Call_RNAseq` | 523 | 239 | 0 | `{"Basal": 78, "Her2": 29, "LumA": 275, "LumB": 124, "Normal": 17, "blank/NA": 239}` |
| `PAM50_mRNA_nature2012` | 207 | 555 | 0 | `{"Basal-like": 35, "HER2-enriched": 14, "Luminal A": 107, "Luminal B": 46, "Normal-like": 5, "blank/NA": 555}` |
| `PR_Status_nature2012` | 457 | 305 | 0 | `{"Indeterminate": 2, "Negative": 140, "Positive": 315, "blank/NA": 305}` |
| `breast_carcinoma_estrogen_receptor_status` | 722 | 40 | 0 | `{"Negative": 162, "Positive": 560, "blank/NA": 40}` |
| `breast_carcinoma_progesterone_receptor_status` | 721 | 41 | 0 | `{"Indeterminate": 2, "Negative": 229, "Positive": 490, "blank/NA": 41}` |
| `her2_immunohistochemistry_level_result` | 392 | 370 | 0 | `{"0": 30, "1+": 172, "2+": 140, "3+": 50, "blank/NA": 370}` |
| `lab_proc_her2_neu_immunohistochemistry_receptor_status` | 627 | 135 | 0 | `{"Equivocal": 139, "Indeterminate": 12, "Negative": 385, "Positive": 91, "blank/NA": 135}` |
| `lab_procedure_her2_neu_in_situ_hybrid_outcome_type` | 258 | 504 | 0 | `{"Equivocal": 4, "Indeterminate": 4, "Negative": 214, "Positive": 36, "blank/NA": 504}` |
| `metastatic_breast_carcinoma_estrogen_receptor_status` | 33 | 729 | 0 | `{"Negative": 12, "Positive": 21, "blank/NA": 729}` |
| `metastatic_breast_carcinoma_progesterone_receptor_status` | 31 | 731 | 0 | `{"Negative": 19, "Positive": 12, "blank/NA": 731}` |

## Matching audit

Rule: `Upper-case sampleID; restrict Xena sample type to 01; match first 12 characters to corrected cohort participant barcode.`  
Participant not found: `0`; blank/NA and nonbinary counts below exclude not-found participants.  
Case mismatches: `0`; participants with multiple primary rows: `0`; maximum primary rows per participant: `1`.

| Phase 14 variable | Not found | Blank/NA | Present nonbinary | Binary available | Sum |
|---|---:|---:|---:|---:|---:|
| ER (`ER_Status_nature2012`) | 0 | 304 | 0 | 458 | 762 |
| PR (`PR_Status_nature2012`) | 0 | 305 | 2 | 455 | 762 |
| HER2 (`HER2_Final_Status_nature2012`) | 0 | 307 | 3 | 452 | 762 |

## Final coalesced coverage

| Receptor | Cohort counts (n=762) | Event counts (n=99) | Conflicts |
|---|---|---|---:|
| ER | `{"Positive": 560, "Negative": 162, "Not available": 40}` | `{"Positive": 61, "Negative": 31, "Not available": 7}` | 0 |
| PR | `{"Positive": 490, "Negative": 229, "Not available": 43}` | `{"Positive": 54, "Negative": 37, "Not available": 8}` | 0 |
| HER2 | `{"Positive": 80, "Negative": 479, "Not available": 203}` | `{"Positive": 10, "Negative": 65, "Not available": 24}` | 163 |

Complete cases (ER, PR, HER2 all Positive/Negative): `530` participants; `72` events.

## Reproducibility

Command: `C:\moenv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase15_task1_receptor_missingness.py --input-dir . --output-dir .\brca-multiomics-survival\results\phase15`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3", "numpy": "2.5.2"}`  
Seeds: `null`.  
Runtime seconds: `4.219178200000897`.
