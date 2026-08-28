#!/usr/bin/env python3
"""Describe available clinical variables for mixed and primary cohorts."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from phase1_common import load_cohort
from prepare_dataset import barcode

FIELDS = {
    "AJCC stage": "ajcc_pathologic_stage.diagnoses",
    "Gender": "gender.demographic",
    "Race": "race.demographic",
    "Ethnicity": "ethnicity.demographic",
}

def summarize(root: Path, cohort: str) -> dict:
    _, y = load_cohort(root, cohort)
    raw = pd.read_csv(root / "TCGA-BRCA.clinical.tsv.gz", sep="\t", low_memory=False)
    raw["sample_key"] = raw["sample"].map(barcode)
    raw = raw.drop_duplicates("sample_key", keep="first").set_index("sample_key").reindex(y.index)
    age = pd.to_numeric(raw["age_at_index.demographic"], errors="coerce")
    out = {"analysis_unit": "record/sample key" if cohort == "mixed" else "unique participant",
           "n": len(raw), "unique_participants": int(pd.Index(raw.index).str[:12].nunique()),
           "age_years": {"nonmissing_n": int(age.notna().sum()), "missing_n": int(age.isna().sum()),
                         "missing_percent": float(age.isna().mean() * 100), "mean": float(age.mean()),
                         "median": float(age.median()), "q1": float(age.quantile(.25)), "q3": float(age.quantile(.75))},
           "categorical": {}}
    for label, col in FIELDS.items():
        values = raw[col]
        out["categorical"][label] = {"source_column": col, "missing_n": int(values.isna().sum()),
                                     "missing_percent": float(values.isna().mean() * 100),
                                     "levels": {str(k): int(v) for k, v in values.dropna().astype(str).value_counts().items()}}
    searched = [c for c in raw.columns if any(x in c.lower() for x in ("estrogen", "progesterone", "her2", "subtype", "pam50"))]
    out["receptor_or_molecular_subtype_columns_found"] = searched
    return out

def main() -> None:
    root = Path(__file__).resolve().parent
    payload = {c: summarize(root, c) for c in ("mixed", "primary")}
    (root / "clinical_characteristics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = ["# Clinical Characteristics", ""]
    for cohort, d in payload.items():
        a = d["age_years"]
        lines += [f"## {cohort}", "", f"Unit: {d['analysis_unit']}; N={d['n']}; unique participants={d['unique_participants']}.", "",
                  f"Age: mean={a['mean']}, median={a['median']}, IQR=[{a['q1']}, {a['q3']}], missing={a['missing_n']} ({a['missing_percent']}%).", ""]
        for label, x in d["categorical"].items():
            lines += [f"### {label}", "", f"Source: `{x['source_column']}`; missing={x['missing_n']} ({x['missing_percent']}%).", "",
                      "| Category | N |", "|---|---:|"] + [f"| {k} | {v} |" for k, v in x["levels"].items()] + [""]
        found = d["receptor_or_molecular_subtype_columns_found"]
        lines += ["ER/PR/HER2 or molecular subtype columns found: " + (", ".join(f"`{x}`" for x in found) if found else "none in the supplied clinical file."), ""]
    (root / "clinical_characteristics.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(payload, indent=2))

if __name__ == "__main__":
    main()
