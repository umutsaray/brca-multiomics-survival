# Supplementary Figure S1 provenance

`figure_s1_calibration_plots.png` was regenerated in Phase 11 by `src/explainability/plot_calibration.py` from the canonical corrected calibration groups in `results/primary_tumor_only/additional_metrics_v2.json`. Those metrics use Phase 7 leakage-corrected DeepSurv OOF risks and unchanged CoxPH/RSF OOF risks.

The six panels cover two cohorts (mixed and corrected primary-tumor-only) and three models (CoxPH, RSF, and DeepSurv). Each panel displays Clinical-only and Full multi-omics calibration series. Points are the stored mean predicted five-year survival and stored Kaplan–Meier five-year survival for each retained risk group; error bars are the stored Greenwood 95% confidence intervals. Point labels report the stored number of unique participants at risk at five years. The diagonal is a visual reference and no new model statistic was computed.

All plotted five-year calibration combinations contain finite retained estimates. Mixed-cohort Full CoxPH contains four rather than five groups because `qcut(..., duplicates="drop")` retained four distinct bins; the figure displays all four without imputation or fabrication. Undefined mixed-cohort one-year AUC values are unrelated to this five-year figure and remain marked unavailable in `results/primary_tumor_only/additional_metrics_v2.md`.
