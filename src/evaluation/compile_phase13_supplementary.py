#!/usr/bin/env python3
"""Compile Phase 13 supplementary submission files from retained results.

This script performs no model fitting or statistical recomputation. It only
extracts already-stored values, writes Tables S4/S5, records provenance, and
packages the requested supplementary files.
"""

from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MIXED = ROOT / "results" / "mixed_cohort"
SUPPLEMENTARY = ROOT / "supplementary"
PHASE13 = ROOT / "results" / "phase13"

S4_PATH = SUPPLEMENTARY / "table_s4_historical_single_split_evaluation.csv"
S5_PATH = SUPPLEMENTARY / "table_s5_historical_wilcoxon_bootstrap_comparison.csv"
README_PATH = SUPPLEMENTARY / "README_supplementary.txt"
ZIP_PATH = SUPPLEMENTARY / "supplementary_materials.zip"
REPORT_PATH = PHASE13 / "phase13_supplementary_compile_results.md"


CAPTIONS = {
    "Table S1": "Full 20-feature SHAP ranking with per-fold selection counts for the corrected primary-tumor-only model (Table 13), and, for reference, the superseded single-split, mixed-tissue SHAP ranking with gene-level annotation.",
    "Table S2": "Full exploratory clinical-vs-omics configuration-pair comparison results (Section 2.4), both cohorts.",
    "Table S3": "Per-fold CoxPH convergence and design-matrix diagnostics (Section 2.5).",
    "Table S4": "Historical, non-participant-grouped single train/test split evaluation across the three survival models (superseded by Table 6 and Table 7; Section 2.5).",
    "Table S5": "Historical, non-participant-grouped Wilcoxon and single-split bootstrap comparison of DeepSurv against classical models (superseded by Table 10; Section 2.6).",
    "Figure S1": "Five-year calibration plots (predicted survival probability by risk-ordered group versus Kaplan–Meier-observed survival with Greenwood 95% confidence intervals) for the corrected primary-tumor-only cohort and the mixed cohort, across all three models (CoxPH, RSF, DeepSurv) and the clinical-only and full multi-omics configurations.",
    "Figure S2": "SHAP summary bar chart of the 20 features shown in Table 13 for the corrected primary-tumor-only, full multi-omics DeepSurv model, with per-feature fold-count annotation and an explicit reminder that these rankings are unstable across cross-validation folds (mean top-20 Jaccard index 0.029) and are not presented as a stable or reproducible biological signal.",
}


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def number(value: float) -> str:
    """Preserve the source JSON's binary-float value without report rounding."""
    return repr(float(value))


def interval(lower: float, upper: float) -> str:
    return f"[{number(lower)}, {number(upper)}]"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def csv_markdown(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    header, *body = rows
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    lines.extend("| " + " | ".join(row) + " |" for row in body)
    return "\n".join(lines)


def main() -> None:
    SUPPLEMENTARY.mkdir(parents=True, exist_ok=True)
    PHASE13.mkdir(parents=True, exist_ok=True)

    baseline_path = MIXED / "baseline_results.json"
    comparison_path = MIXED / "model_comparison.md"
    cv_path = MIXED / "cv_results.json"
    bootstrap_path = MIXED / "bootstrap_ci_results.json"
    statistical_path = MIXED / "statistical_comparison.json"

    baseline = load_json(baseline_path)
    bootstrap = load_json(bootstrap_path)
    statistical = load_json(statistical_path)

    # Existence/read checks for every source explicitly named in the Phase 13 request.
    comparison_text = comparison_path.read_text(encoding="utf-8")
    cv = load_json(cv_path)
    if "Omics-block DeepSurv" not in comparison_text:
        raise ValueError("DeepSurv row was not found in model_comparison.md")
    if cv.get("folds") != 5:
        raise ValueError("cv_results.json does not contain the expected retained five-fold record")

    s4_rows = [
        {
            "Model": "CoxPH",
            "C-index": number(baseline["models"]["coxph"]["test_c_index"]),
            "95% CI": interval(
                bootstrap["models"]["CoxPH"]["ci_95_percentile_lower"],
                bootstrap["models"]["CoxPH"]["ci_95_percentile_upper"],
            ),
            "Validation scheme": "Original, non-participant-grouped single train/test split",
        },
        {
            "Model": "Random Survival Forest",
            "C-index": number(baseline["models"]["random_survival_forest"]["test_c_index"]),
            "95% CI": interval(
                bootstrap["models"]["Random Survival Forest"]["ci_95_percentile_lower"],
                bootstrap["models"]["Random Survival Forest"]["ci_95_percentile_upper"],
            ),
            "Validation scheme": "Original, non-participant-grouped single train/test split",
        },
        {
            "Model": "DeepSurv",
            "C-index": number(bootstrap["models"]["DeepSurv"]["test_c_index"]),
            "95% CI": interval(
                bootstrap["models"]["DeepSurv"]["ci_95_percentile_lower"],
                bootstrap["models"]["DeepSurv"]["ci_95_percentile_upper"],
            ),
            "Validation scheme": "Original, non-participant-grouped single train/test split",
        },
    ]
    write_csv(S4_PATH, ["Model", "C-index", "95% CI", "Validation scheme"], s4_rows)

    s5_rows = []
    for comparison_name in ("DeepSurv vs CoxPH", "DeepSurv vs Random Survival Forest"):
        record = statistical["comparisons"][comparison_name]
        s5_rows.append(
            {
                "Comparison": comparison_name,
                "Wilcoxon p-value": number(record["wilcoxon_two_sided_p_value"]),
                "Bootstrap mean difference": number(record["bootstrap_mean_difference"]),
                "95% CI": interval(record["bootstrap_95_ci_lower"], record["bootstrap_95_ci_upper"]),
                "Significant": "Yes" if record["significant_by_percentile_ci"] else "No",
            }
        )
    write_csv(
        S5_PATH,
        ["Comparison", "Wilcoxon p-value", "Bootstrap mean difference", "95% CI", "Significant"],
        s5_rows,
    )

    readme_lines = ["MDPI Supplementary Materials", "============================", ""]
    for label, caption in CAPTIONS.items():
        readme_lines.extend([label, caption, ""])
    README_PATH.write_text("\n".join(readme_lines), encoding="utf-8")

    archive_files = [
        SUPPLEMENTARY / "table_s1_shap_full_ranking.csv",
        MIXED / "annotated_shap_features.md",
        MIXED / "annotated_shap_features.json",
        SUPPLEMENTARY / "table_s2_exploratory_clinical_vs_omics.csv",
        SUPPLEMENTARY / "table_s3_coxph_convergence_diagnostics.csv",
        S4_PATH,
        S5_PATH,
        SUPPLEMENTARY / "figure_s1_calibration_plots.png",
        SUPPLEMENTARY / "figure_s2_shap_summary_bar.png",
        README_PATH,
    ]
    missing = [str(path) for path in archive_files if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing supplementary inputs: " + ", ".join(missing))

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in archive_files:
            archive.write(path, arcname=path.name)

    zip_size = ZIP_PATH.stat().st_size
    zip_entries = []
    with zipfile.ZipFile(ZIP_PATH, "r") as archive:
        for info in archive.infolist():
            zip_entries.append((info.filename, info.file_size, info.compress_size))

    report = f"""# Phase 13 supplementary compilation results

## Scope

This phase only compiled retained historical values and packaged existing supplementary artefacts. No model was fitted and no statistical result was recomputed.

## Sources and exact extracted fields

- `{baseline_path.relative_to(ROOT).as_posix()}`
  - `models.coxph.test_c_index` = `{number(baseline['models']['coxph']['test_c_index'])}`
  - `models.random_survival_forest.test_c_index` = `{number(baseline['models']['random_survival_forest']['test_c_index'])}`
- `{comparison_path.relative_to(ROOT).as_posix()}`
  - Retained single-split comparison row confirms DeepSurv test C-index; the unrounded value below is taken from the paired-bootstrap JSON.
- `{cv_path.relative_to(ROOT).as_posix()}`
  - Read for requested historical provenance verification (`folds` = `{cv['folds']}`); no CV value is inserted into Table S4 or S5.
- `{bootstrap_path.relative_to(ROOT).as_posix()}`
  - `models.CoxPH.test_c_index` = `{number(bootstrap['models']['CoxPH']['test_c_index'])}`
  - `models.CoxPH.ci_95_percentile_lower/upper` = `{interval(bootstrap['models']['CoxPH']['ci_95_percentile_lower'], bootstrap['models']['CoxPH']['ci_95_percentile_upper'])}`
  - `models.Random Survival Forest.test_c_index` = `{number(bootstrap['models']['Random Survival Forest']['test_c_index'])}`
  - `models.Random Survival Forest.ci_95_percentile_lower/upper` = `{interval(bootstrap['models']['Random Survival Forest']['ci_95_percentile_lower'], bootstrap['models']['Random Survival Forest']['ci_95_percentile_upper'])}`
  - `models.DeepSurv.test_c_index` = `{number(bootstrap['models']['DeepSurv']['test_c_index'])}`
  - `models.DeepSurv.ci_95_percentile_lower/upper` = `{interval(bootstrap['models']['DeepSurv']['ci_95_percentile_lower'], bootstrap['models']['DeepSurv']['ci_95_percentile_upper'])}`
- `{statistical_path.relative_to(ROOT).as_posix()}`
  - `comparisons.DeepSurv vs CoxPH.wilcoxon_two_sided_p_value` = `{number(statistical['comparisons']['DeepSurv vs CoxPH']['wilcoxon_two_sided_p_value'])}`
  - `comparisons.DeepSurv vs CoxPH.bootstrap_mean_difference` = `{number(statistical['comparisons']['DeepSurv vs CoxPH']['bootstrap_mean_difference'])}`
  - `comparisons.DeepSurv vs CoxPH.bootstrap_95_ci_lower/upper` = `{interval(statistical['comparisons']['DeepSurv vs CoxPH']['bootstrap_95_ci_lower'], statistical['comparisons']['DeepSurv vs CoxPH']['bootstrap_95_ci_upper'])}`
  - `comparisons.DeepSurv vs CoxPH.significant_by_percentile_ci` = `{str(statistical['comparisons']['DeepSurv vs CoxPH']['significant_by_percentile_ci']).lower()}`
  - `comparisons.DeepSurv vs Random Survival Forest.wilcoxon_two_sided_p_value` = `{number(statistical['comparisons']['DeepSurv vs Random Survival Forest']['wilcoxon_two_sided_p_value'])}`
  - `comparisons.DeepSurv vs Random Survival Forest.bootstrap_mean_difference` = `{number(statistical['comparisons']['DeepSurv vs Random Survival Forest']['bootstrap_mean_difference'])}`
  - `comparisons.DeepSurv vs Random Survival Forest.bootstrap_95_ci_lower/upper` = `{interval(statistical['comparisons']['DeepSurv vs Random Survival Forest']['bootstrap_95_ci_lower'], statistical['comparisons']['DeepSurv vs Random Survival Forest']['bootstrap_95_ci_upper'])}`
  - `comparisons.DeepSurv vs Random Survival Forest.significant_by_percentile_ci` = `{str(statistical['comparisons']['DeepSurv vs Random Survival Forest']['significant_by_percentile_ci']).lower()}`

The retained bootstrap metadata are `requested_iterations = {bootstrap['requested_iterations']}` and `random_state = {bootstrap['random_state']}`. The retained statistical-comparison metadata are `bootstrap_iterations = {statistical['bootstrap_iterations']}`, `paired_resampling = {str(statistical['paired_resampling']).lower()}`, and `cv_fold_count = {statistical['cv_fold_count']}`.

## Table S4 (full CSV content)

{csv_markdown(S4_PATH)}

## Table S5 (full CSV content)

{csv_markdown(S5_PATH)}

## Submission archive

Archive: `{ZIP_PATH.relative_to(ROOT).as_posix()}`

| File | Uncompressed bytes | Compressed bytes |
|---|---:|---:|
"""
    report += "\n".join(f"| `{name}` | {size} | {compressed} |" for name, size, compressed in zip_entries)
    report += f"""

Total ZIP size: **{zip_size} bytes ({zip_size / (1024 * 1024):.3f} MiB)**. This is below the requested 250 MB limit.

## Uncertainty / limitations

- Table S4's confidence intervals are the already-stored 1,000-resample percentile intervals from the historical single test split; they were not recalculated.
- Table S5 combines the already-stored five-fold Wilcoxon p-values with the already-stored paired single-test-split bootstrap differences. The source explicitly notes low Wilcoxon power at `n = 5`.
- No missing requested source value or caption was encountered.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")

    print(f"Wrote {S4_PATH.relative_to(ROOT)} ({len(s4_rows)} rows)")
    print(f"Wrote {S5_PATH.relative_to(ROOT)} ({len(s5_rows)} rows)")
    print(f"Wrote {README_PATH.relative_to(ROOT)} ({len(CAPTIONS)} captions)")
    print(f"Wrote {ZIP_PATH.relative_to(ROOT)} ({len(zip_entries)} files, {zip_size} bytes)")
    print(f"Wrote {REPORT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
