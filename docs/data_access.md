# Data access

Raw TCGA files are intentionally excluded. The data source is the [UCSC Xena GDC Hub](https://xenabrowser.net/datapages/?hub=https%3A%2F%2Fgdc.xenahubs.net%3A443), hub address `https://gdc.xenahubs.net`. The analysis used these GDC Xena dataset filenames:

- `TCGA-BRCA.star_tpm.tsv.gz`
- `TCGA-BRCA.mirna.tsv.gz`
- `TCGA-BRCA.methylation450.tsv.gz`
- `TCGA-BRCA.clinical.tsv.gz`
- `TCGA-BRCA.survival.tsv.gz`

Direct download pattern: `https://gdc-hub.s3.us-east-1.amazonaws.com/download/<filename>`. For example, methylation is `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.methylation450.tsv.gz`. Dataset availability and releases can change; verify the filenames on the linked hub page before downloading.

Create `data/raw/`, download the five files there, and either run scripts with `--input-dir data/raw` where supported or copy/symlink them into the working input directory expected by the current scripts. The checked-in scripts preserve the exact Phase 1/2 implementation; path-interface harmonization remains a documented TODO. TCGA/GDC data-use terms remain applicable.
