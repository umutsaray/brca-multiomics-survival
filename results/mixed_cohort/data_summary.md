# TCGA-BRCA Multi-Omics Veri Özeti

## Girdiler ve örnek sayıları

- RNA-seq: 1215 benzersiz örnek, 1226 ham örnek sütunu, 60660 özellik
- miRNA: 1189 benzersiz örnek, 1202 ham örnek sütunu, 1881 özellik
- Methylation450: 886 benzersiz örnek, 893 ham örnek sütunu, 486427 özellik
- Clinical: 1242 benzersiz örnek/hasta anahtarı
- Survival: 1219 benzersiz örnek/hasta anahtarı

## Birleştirme

- Omik kesişimindeki nihai örnek sayısı: 853
- Tam tablonun shape'i: 853 x 548975
- Tam çıktı: `merged_multiomics_full.parquet`
- Varyans seçimi: RNA ve methylation blokları için ayrı ayrı en fazla 2000 özellik
- İndirgenmiş tablonun shape'i: 853 x 5888
- İndirgenmiş çıktı: `merged_multiomics_reduced.parquet`

## Eksik veri oranları

- rna: 0.000000%
- mirna: 0.000000%
- meth: 15.539901%

## Sağkalım etiket denetimi

- Saptanan kolonlar: `{'survival_event': 'OS', 'survival_time': 'OS.time', 'clinical_event': 'vital_status.demographic', 'clinical_death_time': 'days_to_death.demographic', 'clinical_followup_time': 'days_to_last_follow_up.diagnoses'}`
- En az bir kaynakta event/time eksiği olan nihai örnek: 11
- Event çelişkisi olan nihai örnek: 0
- Zaman çelişkisi olan nihai örnek: 0
- Clinical yinelenen-grup sütun çatışmaları: 64
- Survival yinelenen-grup sütun çatışmaları: 13
- Eksik veya çelişkili satırlar otomatik olarak silinmedi; bayrak kolonları tam tabloda tutuldu.
