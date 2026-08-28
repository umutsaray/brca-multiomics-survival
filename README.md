# TCGA-BRCA Multi-omics Survival Modeling

This repository contains the reproducible analysis code and derived results for TCGA-BRCA survival prediction using RNA-seq, miRNA, DNA methylation, and clinical covariates. The principal corrected analysis is participant-level and primary-tumor-only, and directly compares clinical-only models with full multi-omics fusion without overstating negative or non-significant findings.

## Repository layout

- `src/data_prep`: dataset assembly and primary-cohort construction.
- `src/preprocessing`: leakage-safe imputation and splitting.
- `src/feature_selection`: train/fold-local univariate Cox selection.
- `src/models`: CoxPH, RSF, DeepSurv, ablation, and Coxnet sensitivity.
- `src/evaluation`: CV, bootstrap, metrics, diagnostics, and experiment logging.
- `src/explainability`: legacy and primary-only SHAP analyses and annotation.
- `results/mixed_cohort`: original mixed-tissue results.
- `results/primary_tumor_only`: corrected primary-tumor-only results.
- `supplementary`: derived supplementary CSV tables and available figures.
- `docs`: data access and manuscript-to-code traceability.
- `archive/superseded`: transparent copies of explicitly superseded artifacts.

Raw TCGA matrices, large design matrices, model binaries, caches, and virtual environments are excluded. No fabricated/example patient data are included.

## Reproduction

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# POSIX: source .venv/bin/activate
pip install -r requirements.txt
```

Download inputs as described in [`docs/data_access.md`](docs/data_access.md). Because the preserved Phase 1/2 scripts import helpers across category directories, set `PYTHONPATH` to all `src` subdirectories or run from an environment that adds them. The logical order is:

1. `prepare_dataset.py`, then `build_primary_cohort.py`
2. `impute_and_split.py`
3. `select_features_v2.py`
4. `ablation_study.py` and `tune_coxph_penalizer.py`
5. evaluation scripts
6. `explain_model_v2.py`

Global split seed is 42. Model-specific seeds and fold participant hashes are recorded row-by-row in [`experiment_log.csv`](experiment_log.csv).

## Traceability

See [`docs/manuscript_table_mapping.md`](docs/manuscript_table_mapping.md). The manuscript file itself was unavailable during assembly, so labels not explicitly supplied by the author are marked for author verification.

## Known limitations

1. No independent METABRIC/GEO external validation cohort was analyzed.
2. The study is retrospective, single-cohort, and has relatively few observed events.
3. Legacy mixed-tissue analyses contain repeated participants and non-primary samples; primary-only results are the principal corrected analysis.
4. Fold-local FDR selection is unstable and can select zero methylation features in some primary folds.
5. Legacy CoxPH runs did not retain Newton iteration counts or inner-fit-level telemetry; convergence warnings persisted.
6. The mixed-cohort DeepSurv internal validation split was not participant-grouped, and preprocessing/selection was not fully repeated inside nested Cox tuning.
7. Stored OOF predictions lack training-fold risk means, preventing retrospective leakage-safe fold centering.
8. ER/PR/HER2 or molecular subtype fields were absent from the supplied clinical file, and fold-level SHAP top-20 stability was low.

The mixed-cohort one-year AUC includes undefined fold values, and the available supplementary calibration/SHAP PNG files come from superseded legacy analyses; both are labelled accordingly.

## License and citation

Code is released under the MIT License. The corresponding-author name remains a placeholder in `LICENSE` because author metadata was not supplied. Citation metadata is TODO pending the final article title, author list, journal citation, and DOI.

## Contact

**TODO:** `[CORRESPONDING AUTHOR NAME] <[CORRESPONDING AUTHOR EMAIL]>`
