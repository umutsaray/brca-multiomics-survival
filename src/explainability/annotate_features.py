#!/usr/bin/env python3
"""Annotate SHAP feature identifiers using MyGene, Illumina 450K, and miRBase."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import shutil
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


MYGENE_URL = "https://mygene.info/v3/gene/{ensembl_id}?fields=symbol,name,entrezgene"
ILLUMINA_MANIFEST_URL = (
    "https://webdata.illumina.com/downloads/productfiles/humanmethylation450/"
    "humanmethylation450_15017482_v1-2.csv"
)
MIRBASE_GFF_URL = "https://www.mirbase.org/download/CURRENT/genomes/hsa.gff3"
NOT_FOUND = "eşleşme bulunamadı"


def parse_shap_table(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pattern = re.compile(
        r"^\|\s*(\d+)\s*\|\s*`([^`]+)`\s*\|\s*([^|]+?)\s*\|\s*([0-9.eE+-]+)\s*\|$"
    )
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            rows.append(
                {
                    "rank": int(match.group(1)),
                    "feature": match.group(2),
                    "block": match.group(3).strip(),
                    "mean_abs_shap": float(match.group(4)),
                }
            )
    if not rows:
        raise ValueError(f"SHAP özellik tablosu ayrıştırılamadı: {path}")
    return rows


def fetch_json(url: str, timeout: int = 60) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "TCGA-BRCA-feature-annotator/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def annotate_rna(ensembl_ids: list[str], cache_path: Path) -> dict[str, dict[str, str]]:
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.is_file() else {}
    result: dict[str, dict[str, str]] = {}
    for raw_id in ensembl_ids:
        query_id = raw_id.split(".", 1)[0]
        if query_id not in cache:
            try:
                payload = fetch_json(MYGENE_URL.format(ensembl_id=urllib.parse.quote(query_id)))
                cache[query_id] = {
                    "symbol": str(payload.get("symbol") or NOT_FOUND),
                    "name": str(payload.get("name") or NOT_FOUND),
                    "source_url": MYGENE_URL.format(ensembl_id=query_id),
                }
            except Exception as exc:
                print(f"UYARI: MyGene sorgusu başarısız ({query_id}): {exc}")
                cache[query_id] = {
                    "symbol": NOT_FOUND,
                    "name": NOT_FOUND,
                    "source_url": MYGENE_URL.format(ensembl_id=query_id),
                }
        result[raw_id] = cache[query_id]
        print(f"MyGene: {raw_id} -> {result[raw_id]['symbol']} | {result[raw_id]['name']}")
    cache_path.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


def download_file(url: str, path: Path) -> None:
    print(f"İndiriliyor: {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "TCGA-BRCA-feature-annotator/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    print(f"İndirme tamamlandı: {path} ({path.stat().st_size} bayt)")


def open_manifest_text(path: Path):
    signature = path.read_bytes()[:4]
    if signature.startswith(b"PK"):
        archive = zipfile.ZipFile(path)
        member = next(name for name in archive.namelist() if name.lower().endswith(".csv"))
        return io.TextIOWrapper(archive.open(member), encoding="utf-8-sig", errors="replace")
    if signature.startswith(b"\x1f\x8b"):
        return gzip.open(path, mode="rt", encoding="utf-8-sig", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def annotate_methylation(
    probes: list[str], manifest_path: Path, manifest_url: str
) -> dict[str, dict[str, str]]:
    if not manifest_path.is_file():
        download_file(manifest_url, manifest_path)
    wanted = set(probes)
    found: dict[str, dict[str, str]] = {}
    with open_manifest_text(manifest_path) as handle:
        header: list[str] | None = None
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            if header is None:
                normalized = [x.strip().lstrip("\ufeff") for x in row]
                if "IlmnID" in normalized and "UCSC_RefGene_Name" in normalized:
                    header = normalized
                    id_index = header.index("IlmnID")
                    gene_index = header.index("UCSC_RefGene_Name")
                continue
            if len(row) <= max(id_index, gene_index):
                continue
            probe = row[id_index].strip()
            if probe in wanted:
                genes = row[gene_index].strip()
                unique_genes = list(dict.fromkeys(x for x in genes.split(";") if x))
                gene_text = ";".join(unique_genes) if unique_genes else NOT_FOUND
                found[probe] = {
                    "symbol": gene_text,
                    "name": gene_text,
                    "source_url": manifest_url,
                }
                print(f"Illumina manifest: {probe} -> {gene_text}")
                if len(found) == len(wanted):
                    break
    if header is None:
        raise ValueError("Manifest içinde IlmnID/UCSC_RefGene_Name başlığı bulunamadı")
    for probe in probes:
        found.setdefault(
            probe,
            {"symbol": NOT_FOUND, "name": NOT_FOUND, "source_url": manifest_url},
        )
        if found[probe]["symbol"] == NOT_FOUND:
            print(f"Illumina manifest: {probe} -> {NOT_FOUND}")
    return found


def parse_gff_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for part in value.strip().split(";"):
        if "=" in part:
            key, item = part.split("=", 1)
            attributes[key] = urllib.parse.unquote(item)
    return attributes


def annotate_mirna(
    ids: list[str], gff_path: Path, reference_path: Path | None = None
) -> dict[str, dict[str, str]]:
    entries: dict[str, dict[str, str]] = {}
    if reference_path and reference_path.is_file():
        with reference_path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row.get("Name", "")
                if name:
                    entries[name.lower()] = row
    missing_from_reference = {item.lower() for item in ids} - set(entries)
    if missing_from_reference and not gff_path.is_file():
        try:
            download_file(MIRBASE_GFF_URL, gff_path)
        except Exception as exc:
            print(f"UYARI: miRBase GFF indirilemedi; yerel referans denenecek: {exc}")
    if gff_path.is_file():
        opener = gzip.open if gff_path.read_bytes()[:2] == b"\x1f\x8b" else open
        with opener(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                attrs = parse_gff_attributes(fields[8])
                name = attrs.get("Name", "")
                if name:
                    entries[name.lower()] = attrs
    result = {}
    for mirna_id in ids:
        attrs = entries.get(mirna_id.lower())
        if attrs:
            accession = attrs.get("ID", NOT_FOUND)
            result[mirna_id] = {
                "symbol": attrs.get("Symbol", attrs.get("Name", mirna_id)),
                "name": attrs.get(
                    "description",
                    f"Homo sapiens {attrs.get('Name', mirna_id)} precursor miRNA; accession {accession}",
                ),
                "family": attrs.get("family") or attrs.get("Family") or NOT_FOUND,
                "source_url": attrs.get("source_url", "https://www.mirbase.org/"),
            }
        else:
            result[mirna_id] = {
                "symbol": NOT_FOUND,
                "name": NOT_FOUND,
                "family": NOT_FOUND,
                "source_url": "https://www.mirbase.org/",
            }
        print(
            f"miRBase: {mirna_id} -> {result[mirna_id]['symbol']} | "
            f"family={result[mirna_id]['family']}"
        )
    return result


def escape_markdown(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--manifest-url", default=ILLUMINA_MANIFEST_URL)
    parser.add_argument("--mirna-reference", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.input_dir.resolve()
    cache_dir = root / "annotation_cache"
    cache_dir.mkdir(exist_ok=True)
    rows = parse_shap_table(root / "shap_top_features.md")
    print(f"SHAP tablosundan {len(rows)} özellik okundu")
    rna_ids = [str(r["feature"]).removeprefix("rna_") for r in rows if r["block"] == "rna"]
    probes = [str(r["feature"]).removeprefix("meth_") for r in rows if r["block"] == "methylation"]
    mirnas = [str(r["feature"]).removeprefix("mirna_") for r in rows if r["block"] == "mirna"]
    rna_annotations = annotate_rna(rna_ids, cache_dir / "mygene_cache.json")
    manifest_path = args.manifest or cache_dir / "HumanMethylation450_15017482_v1-2.csv"
    meth_annotations = annotate_methylation(probes, manifest_path, args.manifest_url)
    mirna_reference = args.mirna_reference or root / "mirbase_reference.tsv"
    mirna_annotations = annotate_mirna(mirnas, cache_dir / "hsa.gff3", mirna_reference)
    output_rows = []
    for row in rows:
        feature = str(row["feature"])
        block = str(row["block"])
        family = "uygulanamaz"
        if block == "rna":
            annotation = rna_annotations[feature.removeprefix("rna_")]
            source = "MyGene.info"
        elif block == "methylation":
            annotation = meth_annotations[feature.removeprefix("meth_")]
            source = "Illumina HumanMethylation450 v1.2 manifest"
        elif block == "mirna":
            annotation = mirna_annotations[feature.removeprefix("mirna_")]
            family = annotation.get("family", NOT_FOUND)
            source = "miRBase"
        else:
            annotation = {
                "symbol": "uygulanamaz",
                "name": feature.removeprefix("clin_").replace("__", ": ").replace("_", " "),
                "source_url": "shap_top_features.md",
            }
            source = "Girdi clinical özellik adı"
        output_rows.append(
            {
                **row,
                "symbol": annotation.get("symbol", NOT_FOUND),
                "name": annotation.get("name", NOT_FOUND),
                "family": family,
                "source": source,
                "source_url": annotation.get("source_url", ""),
            }
        )
    lines = [
        "# Annotated SHAP Features",
        "",
        "Bu tablo yalnız ID→isim eşlemesidir; biyolojik işlev veya hastalık yorumu içermez.",
        "",
        r"| SHAP rank | Original feature | Block | Symbol / associated gene(s) | Official name / annotation | miRNA family | Mean \|SHAP\| | Source |",
        "|---:|---|---|---|---|---|---:|---|",
    ]
    for row in output_rows:
        source_link = f"[{row['source']}]({row['source_url']})" if row["source_url"] else row["source"]
        lines.append(
            f"| {row['rank']} | `{escape_markdown(row['feature'])}` | {row['block']} | "
            f"{escape_markdown(row['symbol'])} | {escape_markdown(row['name'])} | "
            f"{escape_markdown(row['family'])} | {row['mean_abs_shap']:.10g} | {source_link} |"
        )
    (root / "annotated_shap_features.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (root / "annotated_shap_features.json").write_text(
        json.dumps(output_rows, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"annotated_shap_features.md yazıldı: {len(output_rows)} satır")


if __name__ == "__main__":
    main()
