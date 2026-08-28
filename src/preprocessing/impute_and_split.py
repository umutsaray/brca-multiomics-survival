#!/usr/bin/env python3
"""Create leakage-safe TCGA-BRCA train/test matrices.

Order of operations is deliberate: labels are filtered and samples are split
before methylation variance selection or any imputer/preprocessor is fitted.
"""

from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.compose import ColumnTransformer
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from prepare_dataset import barcode, collapse_metadata, find_column


LABEL_COLUMNS = [
    "sample_key",
    "survival_event",
    "survival_time",
    "label_missing_any_source",
]
CLINICAL_COVARIATES = [
    "age_at_index.demographic",
    "gender.demographic",
    "race.demographic",
    "ethnicity.demographic",
    "ajcc_pathologic_stage.diagnoses",
]


def read_schema(path: Path) -> list[str]:
    return pq.ParquetFile(path).schema_arrow.names


def read_columns(path: Path, columns: list[str]) -> pd.DataFrame:
    return pq.read_table(path, columns=columns).to_pandas()


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    serializable = frame.reset_index()
    table = pa.Table.from_pandas(serializable, preserve_index=False).replace_schema_metadata(None)
    pq.write_table(table, path, compression="snappy")


def load_clinical_covariates(path: Path, sample_keys: pd.Index) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    id_col = find_column(header, ["sample", "sample_id", "submitter_id", "_PATIENT", "patient"])
    if id_col is None:
        raise ValueError(f"Clinical örnek ID kolonu bulunamadı. Kolonlar: {header}")
    missing_columns = [c for c in CLINICAL_COVARIATES if c not in header]
    if missing_columns:
        raise ValueError(
            f"Beklenen clinical kovaryat kolonları dosyada yok: {missing_columns}. "
            f"Gerçek kolonlar: {header}"
        )
    print(f"Kullanılan gerçek clinical kovaryatlar: {CLINICAL_COVARIATES}")
    clinical = pd.read_csv(
        path,
        sep="\t",
        usecols=[id_col, *CLINICAL_COVARIATES],
        low_memory=False,
    )
    clinical.index = pd.Index(clinical[id_col].map(barcode), name="sample_key")
    clinical = clinical.drop(columns=id_col)
    clinical, _ = collapse_metadata(clinical, "Clinical kovaryatlar")
    return clinical.reindex(sample_keys)


def scan_train_methylation_variance(
    path: Path,
    train_keys: pd.Index,
    eligible_keys: pd.Index,
    top_n: int,
    chunksize: int,
) -> tuple[list[str], int, float]:
    """Learn the methylation variance ranking using train samples only."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    feature_col, raw_samples = header[0], header[1:]
    selected_raw = [c for c in raw_samples if barcode(c) in eligible_keys]
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=[feature_col, *selected_raw],
        chunksize=chunksize,
        low_memory=False,
    )
    variance_parts: list[pd.Series] = []
    total_probes = 0
    started = time.perf_counter()
    for chunk in reader:
        probe_ids = chunk.pop(feature_col).astype(str)
        values = chunk.apply(pd.to_numeric, errors="coerce").astype("float32")
        values.index = probe_ids
        samples = values.T
        samples.index = pd.Index([barcode(x) for x in samples.index], name="sample_key")
        if samples.index.has_duplicates:
            samples = samples.groupby(level=0, sort=False).mean()
        train_block = samples.reindex(train_keys)
        variances = train_block.var(axis=0, skipna=True, ddof=1)
        variances.index = pd.Index([f"meth_{x}" for x in variances.index])
        variance_parts.append(variances)
        total_probes += len(probe_ids)
        del chunk, values, samples, train_block, variances
        gc.collect()
    ranking = pd.concat(variance_parts).sort_values(ascending=False, na_position="last")
    usable = int(ranking.notna().sum())
    if usable == 0:
        raise ValueError("Train setinde varyansı hesaplanabilen methylation probe'u bulunamadı")
    chosen_n = min(top_n, usable)
    chosen = ranking.iloc[:chosen_n].index.tolist()
    elapsed = time.perf_counter() - started
    print(
        f"Methylation varyans ön filtresi yalnız train setinde öğrenildi: "
        f"{total_probes} probe tarandı, {chosen_n} seçildi, "
        f"{total_probes - chosen_n} probe elendi ({elapsed:.2f} sn)."
    )
    return chosen, total_probes, elapsed


def fit_transform_clinical(
    clinical: pd.DataFrame,
    train_keys: pd.Index,
    test_keys: pd.Index,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    numeric = ["age_at_index.demographic"]
    categorical = [c for c in CLINICAL_COVARIATES if c not in numeric]
    processor = ColumnTransformer(
        transformers=[
            ("numeric", SimpleImputer(strategy="median"), numeric),
            (
                "categorical",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        (
                            "onehot",
                            OneHotEncoder(
                                handle_unknown="ignore",
                                sparse_output=False,
                                dtype=np.float32,
                            ),
                        ),
                    ]
                ),
                categorical,
            ),
        ],
        verbose_feature_names_out=True,
    )
    train_array = processor.fit_transform(clinical.loc[train_keys])
    test_array = processor.transform(clinical.loc[test_keys])
    names = [f"clin_{x}" for x in processor.get_feature_names_out()]
    train = pd.DataFrame(train_array, index=train_keys, columns=names, dtype="float32")
    test = pd.DataFrame(test_array, index=test_keys, columns=names, dtype="float32")
    return train, test


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--top-methylation", type=int, default=10_000)
    parser.add_argument("--knn-neighbors", type=int, default=5)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--methylation-chunksize", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.top_methylation <= 0 or args.knn_neighbors <= 0 or args.methylation_chunksize <= 0:
        raise ValueError("Probe, komşu ve chunk parametreleri pozitif olmalı")
    root = args.input_dir.resolve()
    merged_path = root / "merged_multiomics_full.parquet"
    methylation_path = root / "TCGA-BRCA.methylation450.tsv.gz"
    clinical_path = root / "TCGA-BRCA.clinical.tsv.gz"
    for path in (merged_path, methylation_path, clinical_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    schema = read_schema(merged_path)
    missing_required = [c for c in LABEL_COLUMNS if c not in schema]
    if missing_required:
        raise ValueError(f"Birleşik Parquet'te gerekli kolonlar yok: {missing_required}")
    rna_cols = [c for c in schema if c.startswith("rna_")]
    mirna_cols = [c for c in schema if c.startswith("mirna_")]
    meth_cols = [c for c in schema if c.startswith("meth_")]
    if not rna_cols or not mirna_cols or not meth_cols:
        raise ValueError("RNA, miRNA veya methylation bloklarından en az biri bulunamadı")

    labels = read_columns(merged_path, LABEL_COLUMNS).set_index("sample_key")
    labels.index = labels.index.astype(str)
    initial_count = len(labels)
    remove_mask = labels["label_missing_any_source"].fillna(True).astype(bool)
    removed_count = int(remove_mask.sum())
    labels = labels.loc[~remove_mask, ["survival_event", "survival_time"]].copy()
    labels = labels.rename(columns={"survival_event": "event", "survival_time": "time"})
    if labels[["event", "time"]].isna().any().any():
        raise ValueError("Kontrol bayrağıyla filtreleme sonrası event/time eksiği kaldı")
    if not labels.index.is_unique:
        raise ValueError("sample_key benzersiz değil; hasta bazlı split güvenli değil")
    labels["event"] = labels["event"].astype("int8")
    labels["time"] = pd.to_numeric(labels["time"], errors="raise")
    print(
        f"Eksik sağkalım kontrol bayrağıyla {removed_count} örnek çıkarıldı: "
        f"{initial_count} -> {len(labels)} örnek kaldı."
    )

    train_keys_raw, test_keys_raw = train_test_split(
        labels.index.to_numpy(),
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=labels["event"].to_numpy(),
    )
    train_keys = pd.Index(train_keys_raw, name="sample_key")
    test_keys = pd.Index(test_keys_raw, name="sample_key")
    if len(train_keys.intersection(test_keys)):
        raise RuntimeError("Train/test örnek anahtarları çakışıyor")

    chosen_meth, total_meth, variance_seconds = scan_train_methylation_variance(
        methylation_path,
        train_keys,
        labels.index,
        args.top_methylation,
        args.methylation_chunksize,
    )
    missing_in_parquet = [c for c in chosen_meth if c not in meth_cols]
    if missing_in_parquet:
        raise ValueError(f"Seçilen probe'lar birleşik Parquet'te yok: {missing_in_parquet[:10]}")

    omics = read_columns(merged_path, ["sample_key", *rna_cols, *mirna_cols, *chosen_meth])
    omics = omics.set_index("sample_key").reindex(labels.index)
    rna_nan = int(omics[rna_cols].isna().sum().sum())
    mirna_nan = int(omics[mirna_cols].isna().sum().sum())
    print(f"RNA-seq NaN kontrolü: {rna_nan}")
    print(f"miRNA NaN kontrolü: {mirna_nan}")
    if rna_nan:
        print("UYARI: RNA-seq bloğunda beklenmeyen NaN değerleri var.")
    if mirna_nan:
        print("UYARI: miRNA bloğunda beklenmeyen NaN değerleri var.")

    meth_train = omics.loc[train_keys, chosen_meth]
    meth_test = omics.loc[test_keys, chosen_meth]
    train_missing_before = int(meth_train.isna().sum().sum())
    test_missing_before = int(meth_test.isna().sum().sum())
    imputer = KNNImputer(n_neighbors=args.knn_neighbors, keep_empty_features=True)
    impute_started = time.perf_counter()
    meth_train_array = imputer.fit_transform(meth_train)
    meth_test_array = imputer.transform(meth_test)
    imputation_seconds = time.perf_counter() - impute_started
    meth_train_imputed = pd.DataFrame(
        meth_train_array, index=train_keys, columns=chosen_meth, dtype="float32"
    )
    meth_test_imputed = pd.DataFrame(
        meth_test_array, index=test_keys, columns=chosen_meth, dtype="float32"
    )
    train_missing_after = int(meth_train_imputed.isna().sum().sum())
    test_missing_after = int(meth_test_imputed.isna().sum().sum())
    print(
        f"KNNImputer yalnız train'de fit edildi (k={args.knn_neighbors}); "
        f"train ve test transform süresi {imputation_seconds:.2f} sn. "
        f"NaN: train {train_missing_before}->{train_missing_after}, "
        f"test {test_missing_before}->{test_missing_after}."
    )

    clinical = load_clinical_covariates(clinical_path, labels.index)
    clinical_train, clinical_test = fit_transform_clinical(clinical, train_keys, test_keys)
    base_cols = [*rna_cols, *mirna_cols]
    x_train = pd.concat(
        [omics.loc[train_keys, base_cols], meth_train_imputed, clinical_train], axis=1
    )
    x_test = pd.concat(
        [omics.loc[test_keys, base_cols], meth_test_imputed, clinical_test], axis=1
    )
    y_train = labels.loc[train_keys, ["event", "time"]]
    y_test = labels.loc[test_keys, ["event", "time"]]
    print(f"Nihai X_train shape: {x_train.shape[0]} x {x_train.shape[1]}")
    print(f"Nihai X_test shape: {x_test.shape[0]} x {x_test.shape[1]}")
    train_event_rate = float(y_train["event"].mean())
    test_event_rate = float(y_test["event"].mean())
    overall_event_rate = float(labels["event"].mean())
    print(f"Genel event oranı: {overall_event_rate:.6%}")
    print(f"Train event oranı: {train_event_rate:.6%}")
    print(f"Test event oranı: {test_event_rate:.6%}")

    outputs = {
        "X_train.parquet": x_train,
        "X_test.parquet": x_test,
        "y_train.parquet": y_train,
        "y_test.parquet": y_test,
    }
    for filename, frame in outputs.items():
        write_parquet(frame, root / filename)

    report = [
        "# Leakage-Safe Imputation ve Train/Test Split Özeti",
        "",
        "## Örnekler ve etiketler",
        "",
        f"- Başlangıç örnek sayısı: {initial_count}",
        f"- Eksik survival kontrol bayrağıyla çıkarılan örnek: {removed_count}",
        f"- Kalan örnek sayısı: {len(labels)}",
        f"- Train örnek sayısı: {len(train_keys)}",
        f"- Test örnek sayısı: {len(test_keys)}",
        f"- Genel event oranı: {overall_event_rate:.6%}",
        f"- Train event oranı: {train_event_rate:.6%}",
        f"- Test event oranı: {test_event_rate:.6%}",
        "",
        "## Özellikler ve leakage kontrolleri",
        "",
        f"- RNA özelliği: {len(rna_cols)}; doğrulanan NaN sayısı: {rna_nan}",
        f"- miRNA özelliği: {len(mirna_cols)}; doğrulanan NaN sayısı: {mirna_nan}",
        f"- Başlangıç methylation probe sayısı: {total_meth}",
        f"- Yalnız train varyansıyla seçilen methylation probe sayısı: {len(chosen_meth)}",
        f"- Varyans ön filtresinde elenen probe sayısı: {total_meth - len(chosen_meth)}",
        f"- Varyans tarama süresi: {variance_seconds:.2f} saniye",
        f"- KNNImputer komşu sayısı: {args.knn_neighbors}",
        f"- KNN fit(train) + transform(train/test) süresi: {imputation_seconds:.2f} saniye",
        f"- Methylation NaN train: {train_missing_before} -> {train_missing_after}",
        f"- Methylation NaN test: {test_missing_before} -> {test_missing_after}",
        f"- Clinical kovaryatlar: {CLINICAL_COVARIATES}",
        "- Clinical median/mode imputation ve one-hot encoder yalnız train setinde fit edildi.",
        "- Train/test ayrımı, varyans filtresi ve tüm imputer fit işlemlerinden önce yapıldı.",
        "",
        "## Nihai matrisler",
        "",
        f"- X_train: {x_train.shape[0]} x {x_train.shape[1]}",
        f"- X_test: {x_test.shape[0]} x {x_test.shape[1]}",
        f"- y_train: {y_train.shape[0]} x {y_train.shape[1]}",
        f"- y_test: {y_test.shape[0]} x {y_test.shape[1]}",
        f"- test_size={args.test_size}, random_state={args.random_state}, stratify=event",
    ]
    (root / "split_summary.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n" + "\n".join(report))


if __name__ == "__main__":
    main()
