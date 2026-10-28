#!/usr/bin/env python3
"""Phase 16 Task 0: HER2 audit and clean export of existing Phase 15 tables."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd


IHC = "lab_proc_her2_neu_immunohistochemistry_receptor_status"
ISH = "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"
FINAL = "HER2_Final_Status_nature2012"
BLANK = "blank/NA"


def load(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def scalar(values: list[str]) -> str:
    if not values:
        return BLANK
    if len(values) == 1:
        return str(values[0])
    return " | ".join(sorted(map(str, values)))


def md_table(headers: list[str], rows: list[list[object]]) -> list[str]:
    return [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
        *("| " + " | ".join(str(x) for x in row) + " |" for row in rows),
    ]


def task0(task1: dict, command: str, runtime: float) -> tuple[dict, str]:
    details = task1["coalescing"]["HER2"]["participant_details"]
    records = []
    conflict_tuples: Counter[tuple[str, str, str]] = Counter()
    conflict_pairs: Counter[tuple[str, str, str]] = Counter()
    for participant, item in details.items():
        sources = item["sources"]
        ihc = scalar(sources[IHC])
        ish = scalar(sources[ISH])
        final = scalar(sources[FINAL])
        records.append({"participant": participant, IHC: ihc, ISH: ish, FINAL: final})
        if item["raw_coalesced_result"] == "Conflict":
            conflict_tuples[(ihc, ish, final)] += 1
            values = {IHC: ihc, ISH: ish, FINAL: final}
            names = list(values)
            for i, left in enumerate(names):
                for right in names[i + 1 :]:
                    if values[left] != BLANK and values[right] != BLANK and values[left] != values[right]:
                        conflict_pairs[(left, values[left], right + "=" + values[right])] += 1

    frame = pd.DataFrame(records)
    categories = sorted(set(frame[IHC]) | set(frame[ISH]), key=lambda x: (x == BLANK, x))
    cross = pd.crosstab(frame[IHC], frame[ISH], dropna=False).reindex(index=categories, columns=categories, fill_value=0)
    tuple_rows = []
    for values, count in sorted(conflict_tuples.items(), key=lambda item: (-item[1], item[0])):
        nonblank = [v for v in values if v != BLANK]
        classification = "genuine category disagreement" if len(set(nonblank)) > 1 else "same normalized category / wording only"
        tuple_rows.append({IHC: values[0], ISH: values[1], FINAL: values[2], "n": count, "classification": classification})

    payload = {
        "task": "Phase 16 Task 0(a): HER2 source cross-tabulation and conflict decomposition",
        "status": "completed",
        "decision_rules": {
            "crosstab_sources": [IHC, ISH],
            "blank_label": BLANK,
            "conflict_sources": [IHC, ISH, FINAL],
            "classification": "Distinct normalized nonblank categories are genuine category disagreements; identical normalized categories are wording-only/same-category. No values were recoded.",
            "model_rerun": False,
        },
        "n_participants": len(frame),
        "crosstab": {row: {col: int(cross.loc[row, col]) for col in categories} for row in categories},
        "crosstab_row_totals": {row: int(cross.loc[row].sum()) for row in categories},
        "crosstab_column_totals": {col: int(cross[col].sum()) for col in categories},
        "conflict_participants": int(sum(conflict_tuples.values())),
        "conflict_value_tuples": tuple_rows,
        "pairwise_disagreements_within_conflict_participants": [
            {"source_1": key[0], "value_1": key[1], "source_2_and_value": key[2], "n": count}
            for key, count in sorted(conflict_pairs.items(), key=lambda item: (-item[1], item[0]))
        ],
        "versions": {"python": platform.python_version(), "pandas": pd.__version__},
        "seeds": None,
        "command": command,
        "runtime_seconds": runtime,
    }

    rows = []
    for row in categories:
        rows.append([f"`{row}`", *[int(cross.loc[row, col]) for col in categories], int(cross.loc[row].sum())])
    rows.append(["Column total", *[int(cross[col].sum()) for col in categories], int(cross.to_numpy().sum())])
    lines = [
        "# Phase 16 Task 0(a) — HER2 source cross-tabulation",
        "",
        "## Pre-specified reporting rules",
        "",
        f"- Rows: `{IHC}`; columns: `{ISH}`.",
        f"- `{BLANK}` is explicit. No values were recoded and no model was rerun.",
        "- Distinct normalized nonblank categories are labelled `genuine category disagreement`; identical normalized categories are `same normalized category / wording only`.",
        "- The Phase 15 conflict rule used three sources, so the 163 participants are decomposed by the full IHC/ISH/nature2012 tuple below.",
        "",
        "## IHC × ISH cross-tabulation (all 762 participants)",
        "",
        *md_table(["IHC \\ ISH", *[f"`{x}`" for x in categories], "Row total"], rows),
        "",
        "## Phase 15 HER2 Conflict participants by full value tuple",
        "",
        *md_table(
            ["IHC", "ISH", "HER2 Final nature2012", "n", "classification"],
            [[f"`{r[IHC]}`", f"`{r[ISH]}`", f"`{r[FINAL]}`", r["n"], r["classification"]] for r in tuple_rows],
        ),
        "",
        f"Tuple count sum: `{sum(r['n'] for r in tuple_rows)}`; Phase 15 Conflict count: `{payload['conflict_participants']}`.",
        "",
        "## Reproducibility",
        "",
        f"Command: `{command}`  ",
        f"Versions: `{json.dumps(payload['versions'])}`  ",
        "Seeds: `null`.  ",
        f"Runtime seconds: `{runtime}`.",
    ]
    return payload, "\n".join(lines) + "\n"


def export_phase15_tables(task2a: dict, task2b: dict, command: str) -> str:
    lines = [
        "# Phase 15 tables for manuscript",
        "",
        "This file is a clean export of existing Phase 15 JSON results. No model or statistic was recomputed.",
    ]
    for title, payload in (("Task 2A — all 762 participants", task2a), ("Task 2B — complete cases", task2b)):
        lines.extend(["", f"## {title}", "", "### C-index", ""])
        rows = []
        for config in payload["configurations"]:
            for model in payload["models"]:
                item = payload["summary"][config][model]
                rows.append([config, model, item["mean"], item["sample_sd"], f"`{json.dumps(item['fold_cindices'])}`"])
        lines.extend(md_table(["Configuration", "Model", "Mean", "Sample SD", "Per-fold C-index"], rows))
        lines.extend(["", "### Fold-stratified paired bootstrap", ""])
        rows = []
        for comparison, models in payload["bootstrap_comparisons"].items():
            for model, item in models.items():
                rows.append([
                    comparison,
                    model,
                    item["mean_difference"],
                    f"`{json.dumps(item['ci_95'])}`",
                    item["raw_two_sided_p"],
                    item["bh_fdr_q_within_analysis_family"],
                ])
        lines.extend(md_table(["Comparison", "Model", "Estimate", "95% CI", "Raw p", "BH-adjusted p"], rows))
    lines.extend(["", "## Provenance", "", f"Export command: `{command}`", ""])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    start = time.perf_counter()
    repo = args.repo_dir.resolve()
    phase15 = repo / "results" / "phase15"
    phase16 = repo / "results" / "phase16"
    phase16.mkdir(parents=True, exist_ok=True)
    command = " ".join(sys.argv)
    task1 = load(phase15 / "task1_receptor_missingness.json")
    task2a = load(phase15 / "task2a_missing_indicator.json")
    task2b = load(phase15 / "task2b_complete_case.json")
    payload, report = task0(task1, command, time.perf_counter() - start)
    (phase16 / "task0_her2_crosstab.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (phase16 / "task0_her2_crosstab.md").write_text(report, encoding="utf-8")
    (phase15 / "phase15_tables_for_manuscript.md").write_text(export_phase15_tables(task2a, task2b, command), encoding="utf-8")
    print(json.dumps({"conflicts": payload["conflict_participants"], "tuple_sum": sum(x["n"] for x in payload["conflict_value_tuples"])}, indent=2))


if __name__ == "__main__":
    main()
