#!/usr/bin/env python3
"""Paired statistical comparisons of DeepSurv against survival baselines."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon


COMPARISONS = {
    "DeepSurv vs CoxPH": ("deepsurv", "coxph", "DeepSurv", "CoxPH"),
    "DeepSurv vs Random Survival Forest": (
        "deepsurv",
        "rsf",
        "DeepSurv",
        "Random Survival Forest",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def empirical_two_sided_p_value(differences: np.ndarray) -> float:
    """Finite-sample corrected paired-bootstrap tail probability around zero."""
    n = len(differences)
    lower = (np.count_nonzero(differences <= 0) + 1) / (n + 1)
    upper = (np.count_nonzero(differences >= 0) + 1) / (n + 1)
    return float(min(1.0, 2.0 * min(lower, upper)))


def main() -> None:
    args = parse_args()
    root = args.input_dir.resolve()
    cv = json.loads((root / "cv_results.json").read_text(encoding="utf-8"))
    bootstrap = json.loads(
        (root / "bootstrap_ci_results.json").read_text(encoding="utf-8")
    )
    raw = bootstrap.get("paired_bootstrap_c_indices")
    if raw is None:
        raise ValueError(
            "bootstrap_ci_results.json paired ham dizileri içermiyor. "
            "Önce bootstrap_ci.py'yi aynı seed ile yeniden çalıştırın."
        )
    bootstrap_name = {
        "deepsurv": "DeepSurv",
        "coxph": "CoxPH",
        "rsf": "Random Survival Forest",
    }
    results = {}
    for title, (dl_key, baseline_key, dl_label, baseline_label) in COMPARISONS.items():
        dl_folds = np.asarray(cv["summary"][dl_key]["fold_c_indices"], dtype=float)
        baseline_folds = np.asarray(
            cv["summary"][baseline_key]["fold_c_indices"], dtype=float
        )
        if len(dl_folds) != len(baseline_folds):
            raise ValueError(f"Fold sayıları eşleşmiyor: {title}")
        wilcoxon_result = wilcoxon(
            dl_folds,
            baseline_folds,
            alternative="two-sided",
            zero_method="wilcox",
            method="auto",
        )
        dl_boot = np.asarray(raw[bootstrap_name[dl_key]], dtype=float)
        baseline_boot = np.asarray(raw[bootstrap_name[baseline_key]], dtype=float)
        if len(dl_boot) != len(baseline_boot):
            raise ValueError(f"Paired bootstrap uzunlukları eşleşmiyor: {title}")
        differences = dl_boot - baseline_boot
        lower, upper = np.percentile(differences, [2.5, 97.5])
        includes_zero = bool(lower <= 0 <= upper)
        empirical_p = empirical_two_sided_p_value(differences)
        significant = not includes_zero
        interpretation = (
            "Anlamlı fark bulunamadı; paired bootstrap %95 CI sıfırı içeriyor."
            if includes_zero
            else (
                f"Paired bootstrap %95 CI sıfırı içermiyor; fark {dl_label} lehine."
                if np.mean(differences) > 0
                else f"Paired bootstrap %95 CI sıfırı içermiyor; fark {baseline_label} lehine."
            )
        )
        results[title] = {
            "fold_count": len(dl_folds),
            "wilcoxon_statistic": float(wilcoxon_result.statistic),
            "wilcoxon_two_sided_p_value": float(wilcoxon_result.pvalue),
            "paired_bootstrap_iterations": len(differences),
            "observed_test_c_index_difference": float(
                bootstrap["models"][dl_label]["test_c_index"]
                - bootstrap["models"][baseline_label]["test_c_index"]
            ),
            "bootstrap_mean_difference": float(differences.mean()),
            "bootstrap_median_difference": float(np.median(differences)),
            "bootstrap_95_ci_lower": float(lower),
            "bootstrap_95_ci_upper": float(upper),
            "bootstrap_ci_includes_zero": includes_zero,
            "bootstrap_empirical_two_sided_p_value": empirical_p,
            "significant_by_percentile_ci": significant,
            "interpretation": interpretation,
        }
        print(
            f"{title}: Wilcoxon W={wilcoxon_result.statistic}, p={wilcoxon_result.pvalue}; "
            f"paired bootstrap mean diff={differences.mean():.6f}, "
            f"95% CI=[{lower:.6f}, {upper:.6f}], empirical p={empirical_p:.6f}; "
            f"{interpretation}"
        )

    payload = {
        "cv_fold_count": len(cv["fold_results"]),
        "bootstrap_iterations": bootstrap["requested_iterations"],
        "paired_resampling": True,
        "alpha": args.alpha,
        "comparisons": results,
        "power_note": (
            "Wilcoxon signed-rank test uses only five paired folds (n=5); "
            "statistical power and attainable exact p-value resolution are low."
        ),
    }
    (root / "statistical_comparison.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    lines = [
        "# Statistical Comparison of Survival Models",
        "",
        "The Wilcoxon signed-rank tests use only five paired CV folds (`n=5`). "
        "Therefore statistical power is low and non-significance must not be interpreted as equivalence.",
        "Bootstrap comparisons use the same 1,000 patient resample indices for all three models.",
        "",
        "| Comparison | Wilcoxon W | Wilcoxon p | Observed test ΔC | Bootstrap mean ΔC | Paired-bootstrap 95% CI | Empirical bootstrap p | Conclusion |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for title, values in results.items():
        lines.append(
            f"| {title} | {values['wilcoxon_statistic']:.6g} | "
            f"{values['wilcoxon_two_sided_p_value']:.6g} | "
            f"{values['observed_test_c_index_difference']:+.6f} | "
            f"{values['bootstrap_mean_difference']:+.6f} | "
            f"[{values['bootstrap_95_ci_lower']:+.6f}, {values['bootstrap_95_ci_upper']:+.6f}] | "
            f"{values['bootstrap_empirical_two_sided_p_value']:.6g} | "
            f"{values['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
        ]
    )
    for title, values in results.items():
        lines.append(f"- **{title}:** {values['interpretation']}")
    (root / "statistical_comparison.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
