# Primary-tumor-only Cohort

Selection rule: Keep sample_type 01 only; sort 15-character sample_key lexicographically and keep the first key per 12-character participant. Vials/aliquots were already collapsed at the 15-character key by prepare_dataset.py using their mean.

- source_rows: 842
- source_participants: 769
- source_primary_rows: 762
- primary_participants: 762
- event_1: 99
- censored_0: 663
- removed_from_842: 80
- difference_from_762_primary_rows: 0
- participants_with_multiple_primary_sample_keys: 0
