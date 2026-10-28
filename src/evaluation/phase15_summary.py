#!/usr/bin/env python3
"""Build the Phase 15 summary strictly from completed Phase 15 JSON outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


MODELS = ("coxph", "rsf", "deepsurv")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def ordering_line(label: str, payload: dict, model: str, reference_sign: int | None) -> tuple[str, int]:
    clinical = payload["summary"]["Clinical-only"][model]["mean"]
    full = payload["summary"]["Full"][model]["mean"]
    delta = clinical - full
    sign = 1 if delta > 0 else (-1 if delta < 0 else 0)
    if reference_sign is None:
        changed = "no (this is the unchanged main-analysis ordering)"
    else:
        changed = "yes" if sign != reference_sign else "no"
    return (
        f"- {label}, {model}: Clinical-only={clinical}; Full={full}; "
        f"Clinical-only - Full={delta}; ordering changed={changed}.",
        sign,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase15-dir", type=Path, default=Path("results/phase15"))
    args = parser.parse_args()

    root = args.phase15_dir
    task1 = load_json(root / "task1_receptor_missingness.json")
    task2a = load_json(root / "task2a_missing_indicator.json")
    task2b = load_json(root / "task2b_complete_case.json")

    final = task1["final_coverage"]
    cc = task1["complete_case"]
    lines = [
        "# Phase 15 summary",
        "",
        "All values below are copied or deterministically subtracted from the Phase 15 JSON outputs; no display rounding was applied.",
        "",
        "## Task 1 — ER/PR/HER2 missingness",
        "",
        f"- Status: {task1['status']}.",
        f"- Matched corrected primary-tumor-only participants: {task1['matching']['cohort_participants'] - task1['matching']['participant_not_found']} / {task1['matching']['cohort_participants']}.",
    ]
    for receptor in ("ER", "PR", "HER2"):
        all_counts = final[receptor]["all_762"]
        event_counts = final[receptor]["events_99"]
        conflicts = task1["coalescing"][receptor]["n_conflicts"]
        lines.append(
            f"- {receptor}: all participants={all_counts}; events={event_counts}; conflicts={conflicts}."
        )
    lines.extend(
        [
            f"- Complete-case cohort: n={cc['n_participants']}; events={cc['n_events']}.",
            f"- PAM50 used as model input: {task2a['pam50_used_as_model_input']}.",
            "",
            "## Task 2A — Missing-indicator analysis (all participants)",
            "",
            f"- Status: {task2a['status']}; n={task2a['cohort_n']}; events={task2a['events']}; runtime_seconds={task2a['runtime_seconds']}.",
        ]
    )
    reproduction = task2a["reproduction_check_against_table7"]
    for configuration in ("Clinical-only", "Full"):
        for model in MODELS:
            check = reproduction[configuration][model]
            lines.append(
                f"- Table 7 reproduction, {configuration}, {model}: exact_match={check['exact_match']}; differences={check['differences']}."
            )
    reference_signs: dict[str, int] = {}
    for model in MODELS:
        line, reference_signs[model] = ordering_line("2A", task2a, model, None)
        lines.append(line)

    lines.extend(
        [
            "",
            "## Task 2B — Complete-case analysis",
            "",
            f"- Status: {task2b['status']}; n={task2b['cohort_n']}; events={task2b['events']}; runtime_seconds={task2b['runtime_seconds']}.",
            f"- Outer-test-fold event counts: {task2b['event_counts_per_outer_test_fold']}; folds with fewer than 5 events: {task2b['folds_with_fewer_than_5_test_events']}.",
        ]
    )
    for model in MODELS:
        line, _ = ordering_line("2B (relative to 2A main ordering)", task2b, model, reference_signs[model])
        lines.append(line)

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `task1_receptor_missingness.md` / `.json`",
            "- `task2a_missing_indicator.md` / `.json`",
            "- `task2b_complete_case.md` / `.json`",
        ]
    )
    output = root / "phase15_summary.md"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
