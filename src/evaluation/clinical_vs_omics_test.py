#!/usr/bin/env python3
"""Paired participant-bootstrap comparison of clinical-only and ablation models.

All configurations within a cohort use the same participant resamples.  The
primary, pre-specified comparison is Clinical-only versus Full, separately for
each estimator.  Remaining configurations are explicitly exploratory and their
two-sided bootstrap tail probabilities are Benjamini-Hochberg corrected.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored
from statsmodels.stats.multitest import multipletests

MODELS = ("coxph", "rsf", "deepsurv")

def cindex(frame: pd.DataFrame, risk: str) -> float:
    return float(concordance_index_censored(frame.event.astype(bool), frame.time, frame[risk])[0])

def paired_test(data: pd.DataFrame, iterations: int, seed: int) -> dict:
    configs = data.configuration.drop_duplicates().tolist()
    wide = data.pivot(index=["sample_key", "fold", "event", "time"], columns="configuration", values=list(MODELS))
    wide.columns = [f"{model}::{config}" for model, config in wide.columns]
    wide = wide.reset_index()
    participants = wide.sample_key.str[:12]
    unique = participants.drop_duplicates().to_numpy()
    comparisons = []
    raw = {}
    for model in MODELS:
        clinical = cindex(wide, f"{model}::Clinical-only")
        for config in configs:
            if config == "Clinical-only":
                continue
            key = f"{model}::{config}"
            pair = wide.dropna(subset=[key, f"{model}::Clinical-only"]).copy()
            pair_participants = pair.sample_key.str[:12]
            pair_unique = pair_participants.drop_duplicates().to_numpy()
            point = cindex(pair, key) - cindex(pair, f"{model}::Clinical-only")
            diffs = []
            rng = np.random.default_rng(seed)  # identical index stream for every comparison
            for _ in range(iterations):
                chosen = rng.choice(pair_unique, size=len(pair_unique), replace=True)
                positions = np.concatenate([np.flatnonzero(pair_participants.to_numpy() == pid) for pid in chosen])
                sample = pair.iloc[positions]
                if sample.event.sum() == 0:
                    continue
                diffs.append(cindex(sample, key) - cindex(sample, f"{model}::Clinical-only"))
            arr = np.asarray(diffs)
            lo, hi = np.percentile(arr, [2.5, 97.5])
            pval = min(1.0, 2.0 * min((np.count_nonzero(arr <= 0) + 1) / (len(arr) + 1),
                                      (np.count_nonzero(arr >= 0) + 1) / (len(arr) + 1)))
            label = f"{model}: {config} - Clinical-only"
            raw[label] = arr.tolist()
            comparisons.append({"model": model, "configuration": config,
                                "comparison_role": "pre-specified primary" if config == "Full" else "exploratory",
                                "delta_c_index": point, "ci_95": [float(lo), float(hi)],
                                "contains_zero": bool(lo <= 0 <= hi), "bootstrap_p_value": float(pval),
                                "n_pair_rows": len(pair), "n_pair_participants": len(pair_unique),
                                "valid_bootstrap_iterations": len(arr)})
    exploratory = [i for i, row in enumerate(comparisons) if row["comparison_role"] == "exploratory"]
    adjusted = multipletests([comparisons[i]["bootstrap_p_value"] for i in exploratory], method="fdr_bh")[1]
    for i, q in zip(exploratory, adjusted):
        comparisons[i]["bh_fdr_p_value"] = float(q)
    for row in comparisons:
        if row["comparison_role"] == "pre-specified primary":
            row["bh_fdr_p_value"] = None
        row["interpretation"] = "anlamlı fark bulunamadı" if row["contains_zero"] else "bootstrap %95 GA sıfırı içermiyor"
    return {"n_rows": len(wide), "n_participants": len(unique), "iterations_requested": iterations,
            "seed": seed, "resampling_unit": "12-character TCGA participant ID",
            "comparison_direction": "configuration minus Clinical-only", "comparisons": comparisons,
            "raw_paired_differences": raw}

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--iterations", type=int, default=1000)
    p.add_argument("--random-state", type=int, default=42)
    a = p.parse_args(); root = a.input_dir.resolve(); results = {}
    for cohort, filename in (("mixed", "ablation_oof_extended.parquet"), ("primary", "ablation_oof_primary.parquet")):
        path = root / filename
        if not path.is_file():
            print(f"UYARI: {filename} yok; {cohort} atlandı.")
            continue
        print(f"{cohort}: paired participant bootstrap başlıyor ({a.iterations} iterasyon).")
        results[cohort] = paired_test(pd.read_parquet(path), a.iterations, a.random_state)
    (root / "clinical_vs_omics_test.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    lines = ["# Clinical-only vs Omics: Paired Participant Bootstrap", "",
             "ΔC-index yönü: belirtilen konfigürasyon eksi Clinical-only. Full karşılaştırmaları önceden belirlenmiş primer testlerdir; diğerleri keşifseldir.", ""]
    for cohort, payload in results.items():
        lines += [f"## {cohort}", "", f"Kayıt: {payload['n_rows']}; katılımcı: {payload['n_participants']}; seed: {payload['seed']}.", "",
                  "| Rol | Model | Konfigürasyon | Δ C-index | %95 GA | Ham p | BH-FDR p | Yorum |", "|---|---|---|---:|---|---:|---:|---|"]
        for r in payload["comparisons"]:
            q = "uygulanmadı" if r["bh_fdr_p_value"] is None else str(r["bh_fdr_p_value"])
            lines.append(f"| {r['comparison_role']} | {r['model']} | {r['configuration']} | {r['delta_c_index']} | [{r['ci_95'][0]}, {r['ci_95'][1]}] | {r['bootstrap_p_value']} | {q} | {r['interpretation']} |")
        lines.append("")
    (root / "clinical_vs_omics_test.md").write_text("\n".join(lines), encoding="utf-8")
    print("clinical_vs_omics_test.md/json yazıldı.")

if __name__ == "__main__":
    main()
