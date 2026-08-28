# Event-count reconciliation

## Reproducibility

- Script: `src/evaluation/reconcile_event_counts.py`
- Inputs: `y_train.parquet`, `y_test.parquet`, `primary_tumor_labels.parquet`
- Random seed / folds: not applicable; this is a deterministic census of all retained records.

## Counts

- Mixed cohort records: 842
- Mixed cohort unique participants: 769
- Record-level events: 132
- Participant-level events: 101
- Duplicate-record overcount: 31
- Participants with inconsistent event labels across records: 0
- Participants without a primary-tumor (01) record: 7
- Events among those participants: 2
- Expected events after excluding them: 99
- Events in retained primary-tumor-only labels: 99
- Expected-versus-actual event difference: 0
- Participant label/time mismatches: 0

## Participants without a primary-tumor record

| Participant ID | Event | Sample type(s) | Sample key(s) |
|---|---:|---|---|
| TCGA-A7-A0CE | 0 | 11 | TCGA-A7-A0CE-11 |
| TCGA-A7-A0DB | 0 | 11 | TCGA-A7-A0DB-11 |
| TCGA-BH-A0AY | 0 | 11 | TCGA-BH-A0AY-11 |
| TCGA-BH-A0BV | 0 | 11 | TCGA-BH-A0BV-11 |
| TCGA-BH-A0DZ | 0 | 11 | TCGA-BH-A0DZ-11 |
| TCGA-BH-A18V | 1 | 06 | TCGA-BH-A18V-06 |
| TCGA-BH-A1FD | 1 | 11 | TCGA-BH-A1FD-11 |

## Consistency findings

No participant had conflicting event labels across mixed-cohort records.

## Root cause

The 132-to-99 reduction is fully reconciled: 31 events were duplicate record-level counts, and 2 event-positive participants had no primary-tumor record. The remaining 762 participants have 99 events, exactly matching the primary-tumor-only cohort. No evidence of participant-ID misalignment, inconsistent event labels, or lexicographic-key misassignment was found.
