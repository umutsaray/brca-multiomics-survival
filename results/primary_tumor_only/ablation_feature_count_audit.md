# Ablation-study feature-count audit

## Reproduction settings

- Instrumented script: `src/models/ablation_study.py`
- Cohorts rerun: `mixed` and `primary`
- Outer split: five-fold `StratifiedGroupKFold`, grouped by 12-character participant ID, event-stratified, shuffled with `random_state=42`
- Selection: fold-local univariate Cox; BH-FDR < 0.05; maximum 200 features per omics block
- Model seed: `42 + fold`, unchanged
- Methylation prefilter: 10,000 probes; KNN neighbors: 5
- Reproduced C-index arrays: exactly equal to the retained mixed and primary `ablation_results*.json` result arrays

## Fold-level methylation counts

| Cohort | Fold 1 | Fold 2 | Fold 3 | Fold 4 | Fold 5 |
|---|---:|---:|---:|---:|---:|
| Mixed | 115 | 96 | 159 | 200 | 200 |
| Primary-tumor-only | 0 | 7 | 0 | 1 | 9 |

The count for a fold is shared by every configuration containing that block because ranking occurs once per outer-training fold. Full configuration/block/configuration detail is retained in:

- `results/mixed_cohort/ablation_feature_counts_extended.csv`
- `results/primary_tumor_only/ablation_feature_counts_primary.csv`

## Zero- and low-feature behavior

No mixed-cohort fold selected fewer than five methylation probes. In the primary cohort, folds 1 and 3 selected zero probes, fold 4 selected one, fold 2 selected seven, and fold 5 selected nine.

There is no feature fallback and the BH-FDR threshold is not relaxed. When a single-block configuration has no selected columns, `ablation_study.py` prints `no FDR-qualified features; skipped` and does not fit or append a C-index for that fold. Consequently, primary methylation-only folds 1 and 3 are absent rather than assigned an artificial 0.50. When methylation is empty in a combined configuration, the non-empty clinical/RNA/miRNA columns remain and the model is fitted without methylation. This is why folds 1 and 3 can still contain nonconstant results for configurations whose label includes methylation; in those folds the methylation-containing configuration is algebraically identical to its corresponding configuration without methylation.

The rerun confirmed this behavior directly. For example, in primary fold 1, `Clinical + RNA-seq + Methylation` exactly matched `Clinical + RNA-seq`, `Full` exactly matched `Clinical + RNA-seq + miRNA`, and `Clinical + Methylation` exactly matched `Clinical-only` for all three models.
