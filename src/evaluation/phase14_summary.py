#!/usr/bin/env python3
"""Compile the Phase 14 status summary from the three immutable check outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    a = load(out / "phaseA_er_pr_her2_availability.json")
    b = load(out / "phaseB_rnaseq_scale_sensitivity.json")
    c = load(out / "phaseC_matched_training_size.json")

    av = a["primary_cohort"]["variables"]
    lines = [
        "# Phase 14 quick diagnostic summary",
        "",
        "## Check A — ER / PR / HER2 availability",
        "",
        f"Status: **{a['status']}**.",
        f"ER: available `{av['ER status']['n_available']}`, missing `{av['ER status']['n_missing']}`, events with value `{av['ER status']['n_events_with_variable_available']}` / `{av['ER status']['n_events_total']}`; categories `{json.dumps(av['ER status']['category_counts'])}`.",
        f"PR: available `{av['PR status']['n_available']}`, missing `{av['PR status']['n_missing']}`, events with value `{av['PR status']['n_events_with_variable_available']}` / `{av['PR status']['n_events_total']}`; categories `{json.dumps(av['PR status']['category_counts'])}`.",
        f"HER2: available `{av['HER2 status']['n_available']}`, missing `{av['HER2 status']['n_missing']}`, events with value `{av['HER2 status']['n_events_with_variable_available']}` / `{av['HER2 status']['n_events_total']}`; categories `{json.dumps(av['HER2 status']['category_counts'])}`.",
        f"PAM50: available `{av['PAM50 subtype']['n_available']}`, missing `{av['PAM50 subtype']['n_missing']}`, events with value `{av['PAM50 subtype']['n_events_with_variable_available']}` / `{av['PAM50 subtype']['n_events_total']}`; categories `{json.dumps(av['PAM50 subtype']['category_counts'])}`.",
        "",
        "## Check B — RNA-seq scale",
        "",
        f"Status: **{b['status']}**.",
        b["conclusion"],
        f"Observed minimum `{b['matrix']['statistics_before_pipeline_preprocessing']['minimum']}`, maximum `{b['matrix']['statistics_before_pipeline_preprocessing']['maximum']}`, median `{b['matrix']['statistics_before_pipeline_preprocessing']['median']}`; percentiles `{json.dumps(b['matrix']['statistics_before_pipeline_preprocessing']['percentiles'])}`.",
        "",
        "## Check C — Matched training size",
        "",
        f"Status: **{c['status']}**.",
        "| Configuration | Model | Full mean ± SD | Matched mean ± SD | Matched−full mean difference | 95% CI | p |",
        "|---|---|---|---|---:|---|---:|",
    ]
    for config in ("Clinical-only", "Full"):
        for model in ("coxph", "rsf"):
            rec = c["results"][config][model]
            full, matched, boot = rec["full_outer_training_table7"], rec["matched_85_percent"], rec["matched_minus_full_bootstrap"]
            lines.append(
                f"| {config} | {model} | `{full['mean']} ± {full['sample_sd']}` | `{matched['mean']} ± {matched['sample_sd']}` | "
                f"{boot['observed_mean_difference']} | `{json.dumps(boot['ci_95'])}` | {boot['two_sided_p']} |"
            )
    lines += ["", "Ordering changes:", ""]
    for model, rec in c["ordering_checks"]["clinical_vs_full"].items():
        lines.append(f"- Clinical-only vs Full, `{model}`: `{str(rec['ordering_changed']).lower()}` (full-training delta `{rec['full_training_clinical_minus_full']}`; matched delta `{rec['matched_training_clinical_minus_full']}`).")
    for config, models in c["ordering_checks"]["deepsurv_vs_classical"].items():
        for model, rec in models.items():
            lines.append(f"- `{config}`, DeepSurv vs `{model}`: `{str(rec['ordering_changed']).lower()}` (before `{rec['before_deepsurv_minus_classical']}`; after `{rec['after_deepsurv_minus_classical']}`).")
    lines += [
        "",
        "## Reproducibility records",
        "",
        f"- Check A command: `{a['command']}`; runtime seconds `{a['runtime_seconds']}`; seeds `null`.",
        f"- Check B command: `{b['command']}`; runtime seconds `{b['runtime_seconds']}`; seeds `null`.",
        f"- Check C command: `{c['command']}`; runtime seconds `{c['runtime_seconds']}`; outer seed `{c['method']['outer_seed']}`; bootstrap seed `{c['method']['bootstrap_seed']}`.",
    ]
    (out / "phase14_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out / "phase14_summary.md")


if __name__ == "__main__":
    main()
