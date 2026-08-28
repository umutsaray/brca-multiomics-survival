#!/usr/bin/env python3
"""Prepare TCGA-BRCA multi-omics matrices for survival modelling.

The script keeps clinical and survival labels side by side, reports missing or
conflicting labels, and never drops such rows automatically.
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


BARCODE_LENGTH = 15
FILES = {
    "rna": "TCGA-BRCA.star_tpm.tsv.gz",
    "mirna": "TCGA-BRCA.mirna.tsv.gz",
    "meth": "TCGA-BRCA.methylation450.tsv.gz",
    "clinical": "TCGA-BRCA.clinical.tsv.gz",
    "survival": "TCGA-BRCA.survival.tsv.gz",
}


def barcode(value: object) -> str:
    """Return the requested 15-character TCGA sample key."""
    return str(value).strip()[:BARCODE_LENGTH]


def unique_in_order(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def preview_table(path: Path, label: str) -> None:
    preview = pd.read_csv(path, sep="\t", nrows=5, low_memory=False)
    print(f"\n{label} ilk 5 satır (ilk 12 sütun):")
    print(preview.iloc[:, :12].to_string(index=False))
    print(f"{label} sütunları ({len(preview.columns)}): {list(preview.columns)}")


def read_header(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()


def omics_layout(path: Path, label: str) -> tuple[str, list[str], pd.Index]:
    columns = read_header(path)
    if len(columns) < 2:
        raise ValueError(f"{label}: özellik ID sütunundan sonra örnek sütunu bulunamadı")
    feature_col, sample_cols = columns[0], columns[1:]
    keys = pd.Index(unique_in_order(barcode(x) for x in sample_cols), name="sample_key")
    print(
        f"{label}: {len(sample_cols)} ham örnek sütunu, "
        f"15 karaktere kırpınca {len(keys)} benzersiz örnek"
    )
    return feature_col, sample_cols, keys


def _collapse_duplicate_samples(
    frame: pd.DataFrame, label: str, *, emit_warning: bool = True
) -> pd.DataFrame:
    duplicates = int(frame.index.duplicated(keep=False).sum())
    if duplicates:
        if emit_warning:
            print(
                f"UYARI: {label}: barkod kırpma sonrası {duplicates} satır yinelenen "
                "anahtarlarda; omik değerleri anahtar bazında ortalanıyor."
            )
        frame = frame.groupby(level=0, sort=False).mean()
    return frame


def load_omics(
    path: Path,
    label: str,
    prefix: str,
    allowed_keys: pd.Index,
    chunksize: int | None = None,
) -> pd.DataFrame:
    feature_col, sample_cols, _ = omics_layout(path, label)
    chosen_samples = [c for c in sample_cols if barcode(c) in allowed_keys]
    usecols = [feature_col, *chosen_samples]
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=chunksize,
        low_memory=False,
    )
    chunks = reader if chunksize else [reader]
    blocks: list[pd.DataFrame] = []
    feature_count = 0
    for chunk_number, chunk in enumerate(chunks):
        ids = chunk.pop(feature_col).astype(str)
        # float32 halves the in-memory footprint of the ~450K-feature block and
        # is sufficient for source expression/beta-value precision.
        numeric = chunk.apply(pd.to_numeric, errors="coerce").astype("float32")
        numeric.index = ids
        transposed = numeric.T
        transposed.index = pd.Index([barcode(x) for x in transposed.index], name="sample_key")
        transposed = _collapse_duplicate_samples(
            transposed, label, emit_warning=(chunk_number == 0)
        )
        transposed.columns = [f"{prefix}{x}" for x in ids]
        blocks.append(transposed)
        feature_count += len(ids)
        del chunk, numeric, transposed
        gc.collect()
    result = pd.concat(blocks, axis=1) if len(blocks) > 1 else blocks[0]
    result = result.reindex(allowed_keys)
    print(f"{label}: {result.shape[0]} filtrelenmiş örnek, {feature_count} özellik")
    return result


def _first_non_null(series: pd.Series):
    values = series.dropna()
    return values.iloc[0] if len(values) else np.nan


def collapse_metadata(frame: pd.DataFrame, label: str) -> tuple[pd.DataFrame, int]:
    duplicate_rows = int(frame.index.duplicated(keep=False).sum())
    conflict_count = 0
    if duplicate_rows:
        for _, group in frame.groupby(level=0, sort=False):
            if len(group) > 1:
                conflict_count += sum(group[c].dropna().astype(str).nunique() > 1 for c in group.columns)
        print(
            f"UYARI: {label}: {duplicate_rows} yinelenen satır; ilk eksik olmayan değer "
            f"alındı. Yinelenen gruplarda {conflict_count} sütun-çatışması saptandı."
        )
        frame = frame.groupby(level=0, sort=False).agg(_first_non_null)
    return frame, conflict_count


def find_column(columns: Iterable[str], candidates: list[str]) -> str | None:
    lookup = {str(c).lower(): str(c) for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def load_metadata(path: Path, label: str) -> tuple[pd.DataFrame, int]:
    preview_table(path, label)
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    id_col = find_column(frame.columns, ["sample", "sample_id", "submitter_id", "_PATIENT", "patient"])
    if id_col is None:
        raise ValueError(f"{label}: örnek/hasta kimlik sütunu bulunamadı; görülen sütunlar: {list(frame.columns)}")
    frame.index = pd.Index(frame[id_col].map(barcode), name="sample_key")
    print(f"{label}: {len(frame)} başlangıç satırı, {frame.index.nunique()} benzersiz 15-karakter anahtar")
    return collapse_metadata(frame, label)


def normalize_event(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.lower()
    mapped = text.map(
        {
            "1": 1.0,
            "1.0": 1.0,
            "dead": 1.0,
            "deceased": 1.0,
            "0": 0.0,
            "0.0": 0.0,
            "alive": 0.0,
            "living": 0.0,
        }
    )
    return mapped.astype("Float64")


def build_labels(clinical: pd.DataFrame, survival: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    surv_event_col = find_column(survival.columns, ["OS", "event", "os_event"])
    surv_time_col = find_column(survival.columns, ["OS.time", "time", "os_time"])
    clin_event_col = find_column(clinical.columns, ["OS", "vital_status.demographic", "vital_status"])
    clin_death_col = find_column(clinical.columns, ["days_to_death.demographic", "days_to_death"])
    clin_follow_col = find_column(
        clinical.columns,
        ["days_to_last_follow_up.diagnoses", "days_to_last_follow_up", "days_to_last_known_alive"],
    )
    detected = {
        "survival_event": surv_event_col,
        "survival_time": surv_time_col,
        "clinical_event": clin_event_col,
        "clinical_death_time": clin_death_col,
        "clinical_followup_time": clin_follow_col,
    }
    print(f"Saptanan sağkalım kolonları: {detected}")
    if surv_event_col is None or surv_time_col is None:
        raise ValueError("survival.tsv içinde doğrulanmış event/time sütun çifti bulunamadı")

    labels = pd.DataFrame(index=survival.index.union(clinical.index, sort=False))
    labels["survival_event"] = normalize_event(survival[surv_event_col]).reindex(labels.index)
    labels["survival_time"] = pd.to_numeric(survival[surv_time_col], errors="coerce").reindex(labels.index)
    if clin_event_col:
        labels["clinical_event"] = normalize_event(clinical[clin_event_col]).reindex(labels.index)
    else:
        labels["clinical_event"] = pd.Series(pd.NA, index=labels.index, dtype="Float64")
    death = (
        pd.to_numeric(clinical[clin_death_col], errors="coerce").reindex(labels.index)
        if clin_death_col
        else pd.Series(np.nan, index=labels.index)
    )
    follow = (
        pd.to_numeric(clinical[clin_follow_col], errors="coerce").reindex(labels.index)
        if clin_follow_col
        else pd.Series(np.nan, index=labels.index)
    )
    labels["clinical_time"] = death.where(labels["clinical_event"].eq(1), follow)
    labels["label_missing_any_source"] = labels[
        ["survival_event", "survival_time", "clinical_event", "clinical_time"]
    ].isna().any(axis=1)
    labels["event_conflict"] = (
        labels["survival_event"].notna()
        & labels["clinical_event"].notna()
        & labels["survival_event"].ne(labels["clinical_event"])
    )
    labels["time_conflict"] = (
        labels["survival_time"].notna()
        & labels["clinical_time"].notna()
        & ~np.isclose(labels["survival_time"], labels["clinical_time"], equal_nan=True)
    )
    return labels, detected


def select_top_variance(frame: pd.DataFrame, prefix: str, n: int) -> pd.DataFrame:
    """Return all non-block columns plus the N highest-variance columns in a block."""
    columns = [c for c in frame.columns if c.startswith(prefix)]
    if not columns:
        return frame.copy()
    variances = frame[columns].var(axis=0, skipna=True).sort_values(ascending=False)
    keep = variances.head(min(n, len(variances))).index
    return frame.drop(columns=columns).join(frame.loc[:, keep])


def missing_percent(frame: pd.DataFrame) -> float:
    if frame.size == 0:
        return float("nan")
    return float(frame.isna().to_numpy().mean() * 100)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq

        # Pandas' JSON schema metadata repeats every one of the ~550K column
        # names and can exceed default Thrift reader limits. Keep sample_key as
        # an explicit column and omit that redundant metadata for portability.
        serializable = frame.reset_index()
        table = pa.Table.from_pandas(serializable, preserve_index=False)
        table = table.replace_schema_metadata(None)
        pq.write_table(table, path, compression="snappy")
    except ImportError as exc:
        raise RuntimeError(
            "Parquet yazmak için pyarrow veya fastparquet gerekli. Örn: pip install pyarrow"
        ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--include-methylation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Methylation dosyası varsa dahil et (kapatmak için --no-include-methylation).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Verilirse reduced çıktıda RNA ve methylation bloklarının her birinden en yüksek varyanslı N özellik seçilir.",
    )
    parser.add_argument("--methylation-chunksize", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.top_n is not None and args.top_n <= 0:
        raise ValueError("--top-n pozitif bir tam sayı olmalı")
    root = args.input_dir.resolve()
    paths = {key: root / name for key, name in FILES.items()}
    for required in ("rna", "mirna", "clinical", "survival"):
        if not paths[required].is_file():
            raise FileNotFoundError(paths[required])

    clinical, clinical_dup_conflicts = load_metadata(paths["clinical"], "Clinical")
    survival, survival_dup_conflicts = load_metadata(paths["survival"], "Survival")
    labels, detected = build_labels(clinical, survival)

    layouts: dict[str, tuple[str, list[str], pd.Index]] = {}
    layouts["rna"] = omics_layout(paths["rna"], "RNA-seq")
    layouts["mirna"] = omics_layout(paths["mirna"], "miRNA")
    use_meth = args.include_methylation and paths["meth"].is_file()
    if args.include_methylation and not paths["meth"].is_file():
        print("UYARI: methylation450 dosyası bulunamadı; methylation adımı atlanıyor.")
    elif not args.include_methylation:
        print("Methylation kullanıcı bayrağıyla devre dışı bırakıldı.")
    if use_meth:
        layouts["meth"] = omics_layout(paths["meth"], "Methylation450")

    common = layouts["rna"][2]
    for key in ("mirna", "meth"):
        if key in layouts:
            common = common.intersection(layouts[key][2], sort=False)
    print(f"Omik index kesişimi: {len(common)} gerçek örnek")
    if common.empty:
        raise ValueError("Omik tablolar arasında ortak örnek bulunamadı")

    frames: dict[str, pd.DataFrame] = {}
    frames["rna"] = load_omics(paths["rna"], "RNA-seq", "rna_", common)
    frames["mirna"] = load_omics(paths["mirna"], "miRNA", "mirna_", common)
    if use_meth:
        frames["meth"] = load_omics(
            paths["meth"], "Methylation450", "meth_", common, args.methylation_chunksize
        )

    merged_omics = pd.concat([frames[k] for k in ("rna", "mirna", "meth") if k in frames], axis=1)
    aligned_labels = labels.reindex(common)
    merged = pd.concat([merged_omics, aligned_labels], axis=1)
    full_path = root / "merged_multiomics_full.parquet"
    write_parquet(merged, full_path)

    reduced_path: Path | None = None
    reduced_shape: tuple[int, int] | None = None
    if args.top_n is not None:
        reduced = merged
        reduced = select_top_variance(reduced, "rna_", args.top_n)
        reduced = select_top_variance(reduced, "meth_", args.top_n)
        reduced_path = root / "merged_multiomics_reduced.parquet"
        write_parquet(reduced, reduced_path)
        reduced_shape = reduced.shape

    missing = {key: missing_percent(value) for key, value in frames.items()}
    final_labels = aligned_labels
    missing_rows = int(final_labels["label_missing_any_source"].sum())
    event_conflicts = int(final_labels["event_conflict"].sum())
    time_conflicts = int(final_labels["time_conflict"].sum())
    report = [
        "# TCGA-BRCA Multi-Omics Veri Özeti",
        "",
        "## Girdiler ve örnek sayıları",
        "",
        f"- RNA-seq: {len(layouts['rna'][2])} benzersiz örnek, {len(layouts['rna'][1])} ham örnek sütunu, {frames['rna'].shape[1]} özellik",
        f"- miRNA: {len(layouts['mirna'][2])} benzersiz örnek, {len(layouts['mirna'][1])} ham örnek sütunu, {frames['mirna'].shape[1]} özellik",
    ]
    if "meth" in layouts:
        report.append(
            f"- Methylation450: {len(layouts['meth'][2])} benzersiz örnek, {len(layouts['meth'][1])} ham örnek sütunu, {frames['meth'].shape[1]} özellik"
        )
    else:
        report.append("- Methylation450: dahil edilmedi")
    report.extend(
        [
            f"- Clinical: {len(clinical)} benzersiz örnek/hasta anahtarı",
            f"- Survival: {len(survival)} benzersiz örnek/hasta anahtarı",
            "",
            "## Birleştirme",
            "",
            f"- Omik kesişimindeki nihai örnek sayısı: {len(common)}",
            f"- Tam tablonun shape'i: {merged.shape[0]} x {merged.shape[1]}",
            f"- Tam çıktı: `{full_path.name}`",
        ]
    )
    if reduced_path:
        report.extend(
            [
                f"- Varyans seçimi: RNA ve methylation blokları için ayrı ayrı en fazla {args.top_n} özellik",
                f"- İndirgenmiş tablonun shape'i: {reduced_shape[0]} x {reduced_shape[1]}",
                f"- İndirgenmiş çıktı: `{reduced_path.name}`",
            ]
        )
    else:
        report.append("- İndirgenmiş çıktı üretilmedi (`--top-n` verilmedi; otomatik seçim yapılmadı).")
    report.extend(["", "## Eksik veri oranları", ""])
    for key, percentage in missing.items():
        report.append(f"- {key}: {percentage:.6f}%")
    report.extend(
        [
            "",
            "## Sağkalım etiket denetimi",
            "",
            f"- Saptanan kolonlar: `{detected}`",
            f"- En az bir kaynakta event/time eksiği olan nihai örnek: {missing_rows}",
            f"- Event çelişkisi olan nihai örnek: {event_conflicts}",
            f"- Zaman çelişkisi olan nihai örnek: {time_conflicts}",
            f"- Clinical yinelenen-grup sütun çatışmaları: {clinical_dup_conflicts}",
            f"- Survival yinelenen-grup sütun çatışmaları: {survival_dup_conflicts}",
            "- Eksik veya çelişkili satırlar otomatik olarak silinmedi; bayrak kolonları tam tabloda tutuldu.",
        ]
    )
    report_path = root / "data_summary.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n" + "\n".join(report))


if __name__ == "__main__":
    main()
