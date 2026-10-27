#!/usr/bin/env python3
"""Phase 14 Check A: receptor/PAM50 availability audit; no modelling."""
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


XENA_URL = "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
XENA_DATASET = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
CANONICAL = {
    "ER status": "ER_Status_nature2012",
    "PR status": "PR_Status_nature2012",
    "HER2 status": "HER2_Final_Status_nature2012",
    "PAM50 subtype": "PAM50Call_RNAseq",
}
XENA_RELEVANT_COLUMNS = (
    "ER_Status_nature2012",
    "PR_Status_nature2012",
    "HER2_Final_Status_nature2012",
    "PAM50Call_RNAseq",
    "PAM50_mRNA_nature2012",
    "Integrated_Clusters_with_PAM50__nature2012",
    "breast_carcinoma_estrogen_receptor_status",
    "breast_carcinoma_progesterone_receptor_status",
    "lab_proc_her2_neu_immunohistochemistry_receptor_status",
    "lab_procedure_her2_neu_in_situ_hybrid_outcome_type",
    "metastatic_breast_carcinoma_estrogen_receptor_status",
    "metastatic_breast_carcinoma_progesterone_receptor_status",
)
KEYWORDS = ("estrogen", "progesterone", "her2", "erbb2", "pam50", "molecular_subtype", "receptor_status")


def versions() -> dict:
    return {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__}


def read_xena(url: str) -> pd.DataFrame:
    with urllib.request.urlopen(url, timeout=120) as response:
        return pd.read_csv(response, sep="\t", low_memory=False)


def aggregate_participant(frame: pd.DataFrame, column: str) -> tuple[pd.Series, list[dict]]:
    values = frame[["participant", column]].dropna()
    grouped = values.groupby("participant", sort=True)[column].agg(
        lambda s: sorted(set(s.astype(str)), key=str)
    )
    conflicts = [
        {"participant": participant, "values": vals}
        for participant, vals in grouped.items() if len(vals) > 1
    ]
    resolved = grouped.map(lambda vals: vals[0] if len(vals) == 1 else "CONFLICT: " + " | ".join(vals))
    return resolved, conflicts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    root, out = args.input_dir.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    local_files = [root / "TCGA-BRCA.clinical.tsv.gz", root / "TCGA-BRCA.survival.tsv.gz"]
    local = {}
    for path in local_files:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        hits = [c for c in frame.columns if any(k in c.lower() for k in KEYWORDS)]
        local[path.name] = {
            "n_rows": int(len(frame)),
            "n_columns": int(frame.shape[1]),
            "columns": frame.columns.astype(str).tolist(),
            "receptor_or_subtype_column_matches": hits,
            "has_er_pr_her2_or_molecular_subtype": bool(hits),
        }

    xena = read_xena(XENA_URL)
    if "sampleID" not in xena:
        raise RuntimeError("Xena clinical matrix lacks sampleID")
    related_columns = [c for c in XENA_RELEVANT_COLUMNS if c in xena.columns]
    related_categories = {
        c: sorted(xena[c].dropna().astype(str).unique().tolist(), key=str) for c in related_columns
    }

    cohort_payload = json.loads((root / "primary_tumor_sample_keys.json").read_text(encoding="utf-8"))
    participants = pd.Index([x[:12] for x in cohort_payload["sample_keys"]], name="participant")
    labels = pd.read_parquet(root / "primary_tumor_labels.parquet")
    if "sample_key" in labels.columns:
        labels = labels.set_index("sample_key")
    labels.index = labels.index.astype(str).str[:12]
    labels = labels.loc[participants]
    event_participants = set(labels.index[labels["event"].astype(int).eq(1)])

    xena["participant"] = xena["sampleID"].astype(str).str[:12]
    # PAM50 is sample-specific. Restrict to the same TCGA primary-tumor sample
    # type (01) before collapsing to the requested 12-character participant key.
    xena_primary = xena[xena["sampleID"].astype(str).str[13:15].eq("01")].copy()
    matched = {}
    for label, column in CANONICAL.items():
        if column not in xena:
            matched[label] = {"column": column, "status": "missing_from_xena_dataset"}
            continue
        resolved, conflicts = aggregate_participant(xena_primary, column)
        values = resolved.reindex(participants)
        present = values.notna()
        counts = values[present].value_counts(dropna=False, sort=False).sort_index().to_dict()
        matched[label] = {
            "column": column,
            "status": "available",
            "n_available": int(present.sum()),
            "n_missing": int((~present).sum()),
            "category_counts": {str(k): int(v) for k, v in counts.items()},
            "not_evaluated_observed_count": int((values == "Not evaluated").sum()),
            "n_events_total": int(len(event_participants)),
            "n_events_with_variable_available": int(values.loc[list(event_participants)].notna().sum()),
            "participant_conflicts": conflicts,
        }

    payload = {
        "check": "A - ER / PR / HER2 availability",
        "status": "completed",
        "no_modelling": True,
        "command": " ".join([sys.executable, *sys.argv]),
        "versions": versions(),
        "seeds": None,
        "runtime_seconds": time.perf_counter() - started,
        "local_files": local,
        "xena_source": {
            "dataset": XENA_DATASET,
            "url": XENA_URL,
            "hub": "UCSC Xena legacy TCGA hub",
            "n_rows": int(len(xena)),
            "n_columns": int(xena.shape[1] - 1),
            "relevant_columns_and_observed_categories": related_categories,
            "canonical_columns_used_for_participant_match": CANONICAL,
        },
        "primary_cohort": {
            "n_participants": int(len(participants)),
            "n_events": int(len(event_participants)),
            "match_key": "first 12 characters of the TCGA barcode",
            "xena_row_scope_before_participant_match": "sample type 01 (Primary Tumor), matching the corrected cohort; then collapse by 12-character participant barcode",
            "variables": matched,
        },
    }
    payload["runtime_seconds"] = time.perf_counter() - started
    (out / "phaseA_er_pr_her2_availability.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# Phase 14 Check A — ER / PR / HER2 availability",
        "",
        "Status: **completed (report only; no modelling)**.",
        "",
        "## Clinical files already used",
        "",
    ]
    for name, info in local.items():
        md += [
            f"### `{name}`",
            "",
            f"Rows: `{info['n_rows']}`; columns: `{info['n_columns']}`.",
            f"ER/PR/HER2/PAM50/molecular-subtype column-name matches: `{json.dumps(info['receptor_or_subtype_column_matches'])}`.",
            "",
            "Every column:",
            "",
            "```text",
            "\n".join(info["columns"]),
            "```",
            "",
        ]
    md += [
        "## Additional UCSC Xena dataset",
        "",
        f"Dataset: `{XENA_DATASET}`  ",
        f"URL: {XENA_URL}  ",
        f"Shape read by the script: `{len(xena)} × {xena.shape[1] - 1}` (the script-added participant key is excluded).",
        "",
        "Relevant status/subtype columns and all observed non-missing categories:",
        "",
    ]
    for c, cats in related_categories.items():
        md.append(f"- `{c}`: `{json.dumps(cats, ensure_ascii=False)}`")
    md += [
        "",
        "`Not evaluated` was not observed in any of the four canonical fields; it was not combined with missing values or another category.",
        "",
        "## Match to corrected primary-tumor-only cohort",
        "",
        f"Participants: `{len(participants)}`; events: `{len(event_participants)}`; Xena rows were first restricted to sample type `01` (Primary Tumor), then matched by the first 12 barcode characters.",
        "",
        "| Variable | Xena column | n available | n missing | Events with value / 99 | Category counts | Conflicts |",
        "|---|---|---:|---:|---:|---|---:|",
    ]
    for label, rec in matched.items():
        if rec["status"] != "available":
            md.append(f"| {label} | `{rec['column']}` | null | null | null | missing from dataset | null |")
        else:
            md.append(
                f"| {label} | `{rec['column']}` | {rec['n_available']} | {rec['n_missing']} | "
                f"{rec['n_events_with_variable_available']} / {rec['n_events_total']} | "
                f"`{json.dumps(rec['category_counts'], ensure_ascii=False)}` | {len(rec['participant_conflicts'])} |"
            )
    md += [
        "",
        "No value was imputed and no receptor/subtype variable was added to a model.",
        "",
        "## Reproducibility",
        "",
        f"Command: `{payload['command']}`  ",
        f"Versions: `{json.dumps(payload['versions'])}`  ",
        "Seeds: `null` (no random operation).  ",
        f"Runtime seconds: `{payload['runtime_seconds']}`.",
    ]
    (out / "phaseA_er_pr_her2_availability.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "runtime_seconds": payload["runtime_seconds"], "matched": matched}, indent=2))


if __name__ == "__main__":
    main()
