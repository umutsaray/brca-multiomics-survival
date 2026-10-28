#!/usr/bin/env python3
"""Phase 15 Task 1: pre-specified ER/PR/HER2 missingness and coalescing audit."""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


XENA_DATASET = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
XENA_URL = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"

# Fixed before any Phase 15 model result is inspected.
SOURCE_ORDER = {
    "ER": ["breast_carcinoma_estrogen_receptor_status", "ER_Status_nature2012"],
    "PR": ["breast_carcinoma_progesterone_receptor_status", "PR_Status_nature2012"],
    "HER2": [
        "lab_proc_her2_neu_immunohistochemistry_receptor_status",
        "lab_procedure_her2_neu_in_situ_hybrid_outcome_type",
        "HER2_Final_Status_nature2012",
    ],
}


def is_candidate(name: str) -> bool:
    lower = name.lower()
    if "pam50" in lower:
        return True
    if name in {"ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012"}:
        return True
    receptor = any(token in lower for token in ("estrogen_receptor", "progesterone_receptor", "her2", "erbb"))
    status_like = any(token in lower for token in ("status", "result", "outcome", "call", "final"))
    return receptor and status_like


def normalize(value) -> str | None:
    if pd.isna(value) or not str(value).strip():
        return None
    raw = str(value).strip()
    key = raw.lower().replace("_", " ").replace("-", " ")
    key = " ".join(key.split())
    if key in {"positive", "pos", "+", "3+ positive"}:
        return "Positive"
    if key in {"negative", "neg", "-"}:
        return "Negative"
    if key in {"indeterminate"}:
        return "Indeterminate"
    if key in {"equivocal"}:
        return "Equivocal"
    if key in {"not evaluated", "not evaluated/not available"}:
        return "Not Evaluated"
    if key in {"performed but not available", "performed but not reported"}:
        return "Performed but Not Available"
    return raw


def value_counts_for_cohort(rows: pd.DataFrame, participants: pd.Index, column: str) -> dict:
    grouped = rows.groupby("participant", sort=True)[column].agg(
        lambda s: sorted({str(v).strip() for v in s if pd.notna(v) and str(v).strip()}, key=str)
    )
    counts: dict[str, int] = {}
    nonmissing = 0
    multiple = 0
    for participant in participants:
        values = grouped.get(participant, [])
        if not values:
            counts["blank/NA"] = counts.get("blank/NA", 0) + 1
        elif len(values) == 1:
            nonmissing += 1
            counts[values[0]] = counts.get(values[0], 0) + 1
        else:
            nonmissing += 1
            multiple += 1
            label = "MULTIPLE: " + " | ".join(values)
            counts[label] = counts.get(label, 0) + 1
    return {
        "n_nonmissing": nonmissing,
        "n_blank_or_na": int(len(participants) - nonmissing),
        "value_counts_including_blank_na": dict(sorted(counts.items())),
        "participants_with_multiple_distinct_values": multiple,
    }


def participant_source_values(rows: pd.DataFrame, participants: pd.Index, column: str) -> dict[str, list[str]]:
    grouped = rows.groupby("participant", sort=True)[column].agg(
        lambda s: sorted({normalize(v) for v in s if normalize(v) is not None}, key=str)
    )
    return {p: grouped.get(p, []) for p in participants}


def coalesce(rows: pd.DataFrame, participants: pd.Index, receptor: str) -> tuple[dict[str, str], dict]:
    columns = SOURCE_ORDER[receptor]
    values = {column: participant_source_values(rows, participants, column) for column in columns}
    final: dict[str, str] = {}
    detail: dict[str, dict] = {}
    for participant in participants:
        by_source = {column: values[column][participant] for column in columns}
        all_nonblank = sorted({v for vals in by_source.values() for v in vals}, key=str)
        if len(all_nonblank) > 1:
            raw_result = "Conflict"
            model_value = "Not available"
            chosen_source = None
        else:
            chosen_source = None
            raw_result = None
            for column in columns:
                if by_source[column]:
                    chosen_source = column
                    raw_result = by_source[column][0] if len(by_source[column]) == 1 else "Conflict"
                    break
            if raw_result in {"Positive", "Negative"}:
                model_value = raw_result
            else:
                model_value = "Not available"
        final[participant] = model_value
        detail[participant] = {
            "sources": by_source,
            "chosen_source": chosen_source,
            "raw_coalesced_result": raw_result,
            "model_category": model_value,
        }
    conflict_count = sum(d["raw_coalesced_result"] == "Conflict" for d in detail.values())
    return final, {"participant_details": detail, "n_conflicts": conflict_count}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    root, out = args.input_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    with urllib.request.urlopen(XENA_URL, timeout=120) as response:
        xena = pd.read_csv(response, sep="\t", low_memory=False)
    candidates = [c for c in xena.columns if is_candidate(c)]

    cohort = json.loads((root / "primary_tumor_sample_keys.json").read_text(encoding="utf-8"))
    participants = pd.Index([k[:12] for k in cohort["sample_keys"]], name="participant")
    labels = pd.read_parquet(root / "primary_tumor_labels.parquet")
    if "sample_key" in labels.columns:
        labels = labels.set_index("sample_key")
    labels.index = labels.index.astype(str).str[:12]
    labels = labels.loc[participants]
    events = set(labels.index[labels.event.astype(int).eq(1)])

    xena["sample_id_original"] = xena["sampleID"].astype(str)
    xena["sample_id_upper"] = xena["sample_id_original"].str.upper()
    xena["participant"] = xena["sample_id_upper"].str[:12]
    xena["sample_type"] = xena["sample_id_upper"].str[13:15]
    primary_rows = xena[xena.sample_type.eq("01")].copy()
    found = set(primary_rows.participant)
    not_found = [p for p in participants if p not in found]

    candidate_audit = {c: value_counts_for_cohort(primary_rows, participants, c) for c in candidates}
    coalesced, coalescing_audit = {}, {}
    for receptor in ("ER", "PR", "HER2"):
        coalesced[receptor], coalescing_audit[receptor] = coalesce(primary_rows, participants, receptor)

    final_counts = {}
    for receptor in ("ER", "PR", "HER2"):
        series = pd.Series(coalesced[receptor]).reindex(participants)
        final_counts[receptor] = {
            "all_762": {k: int((series == k).sum()) for k in ("Positive", "Negative", "Not available")},
            "events_99": {k: int((series.loc[list(events)] == k).sum()) for k in ("Positive", "Negative", "Not available")},
            "n_conflicts": coalescing_audit[receptor]["n_conflicts"],
        }

    complete = pd.Series(True, index=participants)
    for receptor in ("ER", "PR", "HER2"):
        complete &= pd.Series(coalesced[receptor]).reindex(participants).isin(["Positive", "Negative"])
    complete_participants = participants[complete]

    missing_breakdown = {}
    phase14_columns = {
        "ER": "ER_Status_nature2012", "PR": "PR_Status_nature2012", "HER2": "HER2_Final_Status_nature2012"
    }
    for receptor, column in phase14_columns.items():
        values = participant_source_values(primary_rows, participants, column)
        nonbinary = 0
        blank = 0
        for p in participants:
            if p in not_found:
                continue
            vals = values[p]
            if not vals:
                blank += 1
            elif not (len(vals) == 1 and vals[0] in {"Positive", "Negative"}):
                nonbinary += 1
        available_binary = sum(len(values[p]) == 1 and values[p][0] in {"Positive", "Negative"} for p in participants)
        missing_breakdown[receptor] = {
            "phase14_column": column,
            "participant_not_found": len(not_found),
            "blank_or_na_among_found": blank,
            "present_but_nonbinary": nonbinary,
            "binary_available": available_binary,
            "sum_check": len(not_found) + blank + nonbinary + available_binary,
        }

    duplicate_primary = primary_rows.groupby("participant").size()
    payload = {
        "task": "Phase 15 Task 1 - receptor missingness",
        "status": "completed",
        "decision_rules_fixed_before_model_results": {
            "source_order": SOURCE_ORDER,
            "fallback": "Use the first nonblank source in the fixed order only when every earlier source is blank.",
            "conflict": "If any available sources contain different normalized nonblank categories, raw result is Conflict and model category is Not available.",
            "binary": "Only Positive and Negative are model-available. Indeterminate, Equivocal, Not Evaluated, Performed but Not Available, Conflict, blank, and every other category are Not available.",
            "pam50": "Availability only; never used as a Phase 15 model input.",
        },
        "command": " ".join([sys.executable, *sys.argv]),
        "versions": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
        "seeds": None,
        "runtime_seconds": None,
        "xena": {
            "dataset": XENA_DATASET, "url": XENA_URL, "rows": int(len(xena)), "columns": int(xena.shape[1] - 4),
            "candidate_columns": candidate_audit,
        },
        "matching": {
            "cohort_participants": len(participants),
            "rule": "Upper-case sampleID; restrict Xena sample type to 01; match first 12 characters to corrected cohort participant barcode.",
            "participant_not_found": len(not_found),
            "participants_not_found": not_found,
            "case_mismatch_count": int((xena.sample_id_original != xena.sample_id_upper).sum()),
            "primary_rows": int(len(primary_rows)),
            "primary_participants": int(primary_rows.participant.nunique()),
            "participants_with_multiple_primary_rows": int((duplicate_primary > 1).sum()),
            "maximum_primary_rows_per_participant": int(duplicate_primary.max()),
        },
        "phase14_missingness_breakdown": missing_breakdown,
        "coalescing": coalescing_audit,
        "final_coverage": final_counts,
        "complete_case": {
            "n_participants": int(complete.sum()),
            "n_events": int(sum(p in events for p in complete_participants)),
            "participants": complete_participants.tolist(),
        },
        "participant_model_categories": {
            p: {receptor: coalesced[receptor][p] for receptor in ("ER", "PR", "HER2")} for p in participants
        },
    }
    payload["runtime_seconds"] = time.perf_counter() - started
    (out / "task1_receptor_missingness.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# Phase 15 Task 1 — ER/PR/HER2 missingness",
        "",
        "## Pre-specified decision rules (fixed before model results)",
        "",
        f"- Source order: `{json.dumps(SOURCE_ORDER)}`.",
        "- Use the first nonblank source only when every earlier source is blank.",
        "- If available sources contain different normalized nonblank categories, mark `Conflict`; do not resolve it.",
        "- Only `Positive` and `Negative` are model-available. Every other result is `Not available`.",
        "- PAM50 is reported only and is not a model input.",
        "",
        "## Candidate columns and counts in the 762-participant cohort",
        "",
        "| Column | Non-missing | Blank/NA | Multiple distinct values | Full value counts |",
        "|---|---:|---:|---:|---|",
    ]
    for column, audit in candidate_audit.items():
        md.append(f"| `{column}` | {audit['n_nonmissing']} | {audit['n_blank_or_na']} | {audit['participants_with_multiple_distinct_values']} | `{json.dumps(audit['value_counts_including_blank_na'], ensure_ascii=False)}` |")
    md += [
        "", "## Matching audit", "",
        f"Rule: `{payload['matching']['rule']}`  ",
        f"Participant not found: `{len(not_found)}`; blank/NA and nonbinary counts below exclude not-found participants.  ",
        f"Case mismatches: `{payload['matching']['case_mismatch_count']}`; participants with multiple primary rows: `{payload['matching']['participants_with_multiple_primary_rows']}`; maximum primary rows per participant: `{payload['matching']['maximum_primary_rows_per_participant']}`.",
        "", "| Phase 14 variable | Not found | Blank/NA | Present nonbinary | Binary available | Sum |", "|---|---:|---:|---:|---:|---:|",
    ]
    for receptor, rec in missing_breakdown.items():
        md.append(f"| {receptor} (`{rec['phase14_column']}`) | {rec['participant_not_found']} | {rec['blank_or_na_among_found']} | {rec['present_but_nonbinary']} | {rec['binary_available']} | {rec['sum_check']} |")
    md += ["", "## Final coalesced coverage", "", "| Receptor | Cohort counts (n=762) | Event counts (n=99) | Conflicts |", "|---|---|---|---:|"]
    for receptor, rec in final_counts.items():
        md.append(f"| {receptor} | `{json.dumps(rec['all_762'])}` | `{json.dumps(rec['events_99'])}` | {rec['n_conflicts']} |")
    md += [
        "",
        f"Complete cases (ER, PR, HER2 all Positive/Negative): `{payload['complete_case']['n_participants']}` participants; `{payload['complete_case']['n_events']}` events.",
        "",
        "## Reproducibility",
        "",
        f"Command: `{payload['command']}`  ",
        f"Versions: `{json.dumps(payload['versions'])}`  ",
        "Seeds: `null`.  ",
        f"Runtime seconds: `{payload['runtime_seconds']}`.",
    ]
    (out / "task1_receptor_missingness.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "final_coverage": final_counts, "complete_case": payload["complete_case"], "runtime_seconds": payload["runtime_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
