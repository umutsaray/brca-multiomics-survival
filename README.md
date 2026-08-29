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

See [`docs/manuscript_table_mapping.md`](docs/manuscript_table_mapping.md) for the complete final-order mapping of every manuscript table and supplementary figure to its producing code and retained output.

## Known limitations

1. No independent METABRIC/GEO external validation cohort was analyzed.
2. The study is retrospective, single-cohort, and has relatively few observed events.
3. Legacy mixed-tissue analyses contain repeated participants and non-primary samples; primary-only results are the principal corrected analysis.
4. Fold-local FDR selection is unstable and can select zero methylation features in some primary folds.
5. Legacy CoxPH runs did not retain Newton iteration counts or inner-fit-level telemetry; convergence warnings persisted.
6. Phases 7 and 9 corrected DeepSurv's internal-validation leakage across all ten ablation configurations. Five Methylation-only fold attempts had zero internal-train FDR-qualified probes and remain explicit missing estimates rather than using fallback features or historical risks.
7. Raw risks are not comparable across outer folds. Phase 7 therefore uses within-fold concordance and participant-level fold-stratified bootstrap rather than pooled OOF concordance for the primary comparison.
8. ER/PR/HER2 or molecular subtype fields were absent from the supplied clinical file, and fold-level SHAP top-20 stability remained low after leakage correction.

The mixed-cohort one-year AUC includes undefined fold values and remains explicitly reported as unavailable. Supplementary calibration and SHAP figures were regenerated from the Phase 7 corrected DeepSurv outputs.

## License and citation

Code is released under the MIT License. Citation metadata will be added after the final article title, author list, journal citation, and DOI are available.

## Contact

Umut Saray, [umutsaray@gmail.com](mailto:umutsaray@gmail.com)
