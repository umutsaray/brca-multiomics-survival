# Phase 14 Check A — ER / PR / HER2 availability

Status: **completed (report only; no modelling)**.

## Clinical files already used

### `TCGA-BRCA.clinical.tsv.gz`

Rows: `1255`; columns: `85`.
ER/PR/HER2/PAM50/molecular-subtype column-name matches: `[]`.

Every column:

```text
sample
id
disease_type
case_id
submitter_id
primary_site
alcohol_history.exposures
race.demographic
gender.demographic
ethnicity.demographic
vital_status.demographic
age_at_index.demographic
days_to_birth.demographic
year_of_birth.demographic
year_of_death.demographic
primary_site.project
project_id.project
disease_type.project
name.project
name.program.project
tissue_source_site_id.tissue_source_site
code.tissue_source_site
name.tissue_source_site
project.tissue_source_site
bcr_id.tissue_source_site
days_to_death.demographic
entity_submitter_id.annotations
notes.annotations
submitter_id.annotations
classification.annotations
entity_id.annotations
created_datetime.annotations
annotation_id.annotations
entity_type.annotations
updated_datetime.annotations
case_id.annotations
state.annotations
category.annotations
status.annotations
case_submitter_id.annotations
synchronous_malignancy.diagnoses
ajcc_pathologic_stage.diagnoses
days_to_diagnosis.diagnoses
last_known_disease_status.diagnoses
tissue_or_organ_of_origin.diagnoses
days_to_last_follow_up.diagnoses
age_at_diagnosis.diagnoses
primary_diagnosis.diagnoses
prior_malignancy.diagnoses
year_of_diagnosis.diagnoses
prior_treatment.diagnoses
ajcc_staging_system_edition.diagnoses
ajcc_pathologic_t.diagnoses
morphology.diagnoses
ajcc_pathologic_n.diagnoses
ajcc_pathologic_m.diagnoses
classification_of_tumor.diagnoses
icd_10_code.diagnoses
site_of_resection_or_biopsy.diagnoses
tumor_grade.diagnoses
progression_or_recurrence.diagnoses
age_at_earliest_diagnosis.diagnoses.xena_derived
age_at_earliest_diagnosis_in_years.diagnoses.xena_derived
treatment_id.treatments.diagnoses
submitter_id.treatments.diagnoses
treatment_type.treatments.diagnoses
treatment_or_therapy.treatments.diagnoses
created_datetime.treatments.diagnoses
updated_datetime.treatments.diagnoses
state.treatments.diagnoses
sample_type_id.samples
tumor_descriptor.samples
sample_id.samples
sample_type.samples
composition.samples
days_to_collection.samples
initial_weight.samples
preservation_method.samples
pathology_report_uuid.samples
oct_embedded.samples
specimen_type.samples
days_to_sample_procurement.samples
is_ffpe.samples
tissue_type.samples
annotations.samples
```

### `TCGA-BRCA.survival.tsv.gz`

Rows: `1232`; columns: `4`.
ER/PR/HER2/PAM50/molecular-subtype column-name matches: `[]`.

Every column:

```text
sample
OS.time
OS
_PATIENT
```

## Additional UCSC Xena dataset

Dataset: `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`  
URL: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix  
Shape read by the script: `1247 × 194` (the script-added participant key is excluded).

Relevant status/subtype columns and all observed non-missing categories:

- `ER_Status_nature2012`: `["Indeterminate", "Negative", "Positive"]`
- `PR_Status_nature2012`: `["Indeterminate", "Negative", "Positive"]`
- `HER2_Final_Status_nature2012`: `["Equivocal", "Negative", "Positive"]`
- `PAM50Call_RNAseq`: `["Basal", "Her2", "LumA", "LumB", "Normal"]`
- `PAM50_mRNA_nature2012`: `["Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like"]`
- `Integrated_Clusters_with_PAM50__nature2012`: `["1.0", "2.0", "3.0", "4.0"]`
- `breast_carcinoma_estrogen_receptor_status`: `["Indeterminate", "Negative", "Positive"]`
- `breast_carcinoma_progesterone_receptor_status`: `["Indeterminate", "Negative", "Positive"]`
- `lab_proc_her2_neu_immunohistochemistry_receptor_status`: `["Equivocal", "Indeterminate", "Negative", "Positive"]`
- `lab_procedure_her2_neu_in_situ_hybrid_outcome_type`: `["Equivocal", "Indeterminate", "Negative", "Positive"]`
- `metastatic_breast_carcinoma_estrogen_receptor_status`: `["Negative", "Positive"]`
- `metastatic_breast_carcinoma_progesterone_receptor_status`: `["Negative", "Positive"]`

`Not evaluated` was not observed in any of the four canonical fields; it was not combined with missing values or another category.

## Match to corrected primary-tumor-only cohort

Participants: `762`; events: `99`; Xena rows were first restricted to sample type `01` (Primary Tumor), then matched by the first 12 barcode characters.

| Variable | Xena column | n available | n missing | Events with value / 99 | Category counts | Conflicts |
|---|---|---:|---:|---:|---|---:|
| ER status | `ER_Status_nature2012` | 458 | 304 | 64 / 99 | `{"Negative": 103, "Positive": 355}` | 0 |
| PR status | `PR_Status_nature2012` | 457 | 305 | 64 / 99 | `{"Indeterminate": 2, "Negative": 140, "Positive": 315}` | 0 |
| HER2 status | `HER2_Final_Status_nature2012` | 455 | 307 | 59 / 99 | `{"Equivocal": 3, "Negative": 395, "Positive": 57}` | 0 |
| PAM50 subtype | `PAM50Call_RNAseq` | 523 | 239 | 73 / 99 | `{"Basal": 78, "Her2": 29, "LumA": 275, "LumB": 124, "Normal": 17}` | 0 |

No value was imputed and no receptor/subtype variable was added to a model.

## Reproducibility

Command: `C:\Users\UMUT\Desktop\Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer\.venv\Scripts\python.exe .\brca-multiomics-survival\src\evaluation\phase14_a_clinical_availability.py --input-dir . --output-dir .\brca-multiomics-survival\results\phase14`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3", "numpy": "2.5.2"}`  
Seeds: `null` (no random operation).  
Runtime seconds: `2.6271913999225944`.
