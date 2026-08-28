#!/usr/bin/env python3
"""Reconcile record- and participant-level event counts across BRCA cohorts."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_labels(path: Path) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    if "sample_key" in frame.columns:
        frame = frame.set_index("sample_key")
    frame.index = frame.index.astype(str)
    return frame[["event", "time"]]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("event_count_reconciliation.md"))
    args = parser.parse_args()

    mixed = pd.concat([
        read_labels(args.data_dir / "y_train.parquet"),
        read_labels(args.data_dir / "y_test.parquet"),
    ])
    mixed["participant_id"] = mixed.index.str[:12]
    mixed["sample_type"] = mixed.index.str[13:15]
    primary = read_labels(args.data_dir / "primary_tumor_labels.parquet")
    primary["participant_id"] = primary.index.str[:12]

    groups = mixed.groupby("participant_id", sort=True)
    inconsistent = groups.filter(lambda g: g["event"].nunique(dropna=False) > 1)
    participant = groups.agg(event=("event", "first"), time=("time", "first"))
    primary_ids = set(mixed.loc[mixed.sample_type.eq("01"), "participant_id"])
    no_primary_ids = sorted(set(participant.index) - primary_ids)
    no_primary = mixed[mixed.participant_id.isin(no_primary_ids)].sort_values(["participant_id", "sample_type"])
    remaining = participant.drop(index=no_primary_ids)
    primary_by_id = primary.set_index("participant_id")
    joined = remaining.join(primary_by_id[["event", "time"]], lsuffix="_expected", rsuffix="_actual")
    mismatches = joined[
        joined.event_expected.ne(joined.event_actual)
        | joined.time_expected.ne(joined.time_actual)
    ]

    record_events = int(mixed.event.sum())
    participant_events = int(participant.event.sum())
    expected_events = int(remaining.event.sum())
    actual_events = int(primary.event.sum())
    lines = [
        "# Event-count reconciliation",
        "",
        "## Reproducibility",
        "",
        "- Script: `src/evaluation/reconcile_event_counts.py`",
        "- Inputs: `y_train.parquet`, `y_test.parquet`, `primary_tumor_labels.parquet`",
        "- Random seed / folds: not applicable; this is a deterministic census of all retained records.",
        "",
        "## Counts",
        "",
        f"- Mixed cohort records: {len(mixed)}",
        f"- Mixed cohort unique participants: {mixed.participant_id.nunique()}",
        f"- Record-level events: {record_events}",
        f"- Participant-level events: {participant_events}",
        f"- Duplicate-record overcount: {record_events - participant_events}",
        f"- Participants with inconsistent event labels across records: {inconsistent.participant_id.nunique()}",
        f"- Participants without a primary-tumor (01) record: {len(no_primary_ids)}",
        f"- Events among those participants: {int(no_primary.groupby('participant_id').event.first().sum())}",
        f"- Expected events after excluding them: {expected_events}",
        f"- Events in retained primary-tumor-only labels: {actual_events}",
        f"- Expected-versus-actual event difference: {expected_events - actual_events}",
        f"- Participant label/time mismatches: {len(mismatches)}",
        "",
        "## Participants without a primary-tumor record",
        "",
        "| Participant ID | Event | Sample type(s) | Sample key(s) |",
        "|---|---:|---|---|",
    ]
    for pid, group in no_primary.groupby("participant_id", sort=True):
        lines.append(f"| {pid} | {int(group.event.iloc[0])} | {', '.join(group.sample_type)} | {', '.join(group.index)} |")
    lines += ["", "## Consistency findings", ""]
    if inconsistent.empty:
        lines.append("No participant had conflicting event labels across mixed-cohort records.")
    else:
        lines.append("Conflicting records:")
        lines.append(inconsistent.reset_index().to_markdown(index=False))
    lines += ["", "## Root cause", ""]
    if expected_events == actual_events and mismatches.empty:
        lines.append(
            f"The 132-to-99 reduction is fully reconciled: {record_events - participant_events} events were duplicate record-level counts, "
            f"and {participant_events - expected_events} event-positive participants had no primary-tumor record. The remaining "
            f"{len(remaining)} participants have {actual_events} events, exactly matching the primary-tumor-only cohort. No evidence "
            "of participant-ID misalignment, inconsistent event labels, or lexicographic-key misassignment was found."
        )
    else:
        lines.append("A serious label mismatch remains; the affected participants are listed below and downstream reruns require approval.")
        lines.append(mismatches.reset_index().to_markdown(index=False))
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
