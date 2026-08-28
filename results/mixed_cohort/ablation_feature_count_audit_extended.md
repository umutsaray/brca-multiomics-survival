# Ablation feature-count audit

- Cohort: mixed
- Outer folds: 5-fold StratifiedGroupKFold, participant-grouped, event-stratified, shuffle seed 42.
- Selection: fold-local univariate Cox, BH-FDR < 0.05, capped at 200 per omics block.
- No FDR-threshold relaxation or feature fallback exists. A zero-feature single-block configuration is skipped; in a combined configuration, the empty block contributes no columns and the remaining blocks are fitted.

| Configuration | Fold | Block | Selected features |
|---|---:|---|---:|
| Clinical-only | 1 | clinical | 15 |
| RNA-seq only | 1 | rna | 200 |
| Methylation only | 1 | methylation | 115 |
| miRNA only | 1 | mirna | 19 |
| Clinical + RNA-seq | 1 | rna | 200 |
| Clinical + RNA-seq | 1 | clinical | 15 |
| Clinical + RNA-seq + Methylation | 1 | rna | 200 |
| Clinical + RNA-seq + Methylation | 1 | methylation | 115 |
| Clinical + RNA-seq + Methylation | 1 | clinical | 15 |
| Clinical + RNA-seq + miRNA | 1 | rna | 200 |
| Clinical + RNA-seq + miRNA | 1 | mirna | 19 |
| Clinical + RNA-seq + miRNA | 1 | clinical | 15 |
| Full | 1 | rna | 200 |
| Full | 1 | mirna | 19 |
| Full | 1 | methylation | 115 |
| Full | 1 | clinical | 15 |
| Clinical + Methylation | 1 | methylation | 115 |
| Clinical + Methylation | 1 | clinical | 15 |
| Clinical + miRNA | 1 | mirna | 19 |
| Clinical + miRNA | 1 | clinical | 15 |
| Clinical-only | 2 | clinical | 16 |
| RNA-seq only | 2 | rna | 200 |
| Methylation only | 2 | methylation | 96 |
| miRNA only | 2 | mirna | 9 |
| Clinical + RNA-seq | 2 | rna | 200 |
| Clinical + RNA-seq | 2 | clinical | 16 |
| Clinical + RNA-seq + Methylation | 2 | rna | 200 |
| Clinical + RNA-seq + Methylation | 2 | methylation | 96 |
| Clinical + RNA-seq + Methylation | 2 | clinical | 16 |
| Clinical + RNA-seq + miRNA | 2 | rna | 200 |
| Clinical + RNA-seq + miRNA | 2 | mirna | 9 |
| Clinical + RNA-seq + miRNA | 2 | clinical | 16 |
| Full | 2 | rna | 200 |
| Full | 2 | mirna | 9 |
| Full | 2 | methylation | 96 |
| Full | 2 | clinical | 16 |
| Clinical + Methylation | 2 | methylation | 96 |
| Clinical + Methylation | 2 | clinical | 16 |
| Clinical + miRNA | 2 | mirna | 9 |
| Clinical + miRNA | 2 | clinical | 16 |
| Clinical-only | 3 | clinical | 16 |
| RNA-seq only | 3 | rna | 200 |
| Methylation only | 3 | methylation | 159 |
| miRNA only | 3 | mirna | 28 |
| Clinical + RNA-seq | 3 | rna | 200 |
| Clinical + RNA-seq | 3 | clinical | 16 |
| Clinical + RNA-seq + Methylation | 3 | rna | 200 |
| Clinical + RNA-seq + Methylation | 3 | methylation | 159 |
| Clinical + RNA-seq + Methylation | 3 | clinical | 16 |
| Clinical + RNA-seq + miRNA | 3 | rna | 200 |
| Clinical + RNA-seq + miRNA | 3 | mirna | 28 |
| Clinical + RNA-seq + miRNA | 3 | clinical | 16 |
| Full | 3 | rna | 200 |
| Full | 3 | mirna | 28 |
| Full | 3 | methylation | 159 |
| Full | 3 | clinical | 16 |
| Clinical + Methylation | 3 | methylation | 159 |
| Clinical + Methylation | 3 | clinical | 16 |
| Clinical + miRNA | 3 | mirna | 28 |
| Clinical + miRNA | 3 | clinical | 16 |
| Clinical-only | 4 | clinical | 16 |
| RNA-seq only | 4 | rna | 200 |
| Methylation only | 4 | methylation | 200 |
| miRNA only | 4 | mirna | 23 |
| Clinical + RNA-seq | 4 | rna | 200 |
| Clinical + RNA-seq | 4 | clinical | 16 |
| Clinical + RNA-seq + Methylation | 4 | rna | 200 |
| Clinical + RNA-seq + Methylation | 4 | methylation | 200 |
| Clinical + RNA-seq + Methylation | 4 | clinical | 16 |
| Clinical + RNA-seq + miRNA | 4 | rna | 200 |
| Clinical + RNA-seq + miRNA | 4 | mirna | 23 |
| Clinical + RNA-seq + miRNA | 4 | clinical | 16 |
| Full | 4 | rna | 200 |
| Full | 4 | mirna | 23 |
| Full | 4 | methylation | 200 |
| Full | 4 | clinical | 16 |
| Clinical + Methylation | 4 | methylation | 200 |
| Clinical + Methylation | 4 | clinical | 16 |
| Clinical + miRNA | 4 | mirna | 23 |
| Clinical + miRNA | 4 | clinical | 16 |
| Clinical-only | 5 | clinical | 16 |
| RNA-seq only | 5 | rna | 200 |
| Methylation only | 5 | methylation | 200 |
| miRNA only | 5 | mirna | 7 |
| Clinical + RNA-seq | 5 | rna | 200 |
| Clinical + RNA-seq | 5 | clinical | 16 |
| Clinical + RNA-seq + Methylation | 5 | rna | 200 |
| Clinical + RNA-seq + Methylation | 5 | methylation | 200 |
| Clinical + RNA-seq + Methylation | 5 | clinical | 16 |
| Clinical + RNA-seq + miRNA | 5 | rna | 200 |
| Clinical + RNA-seq + miRNA | 5 | mirna | 7 |
| Clinical + RNA-seq + miRNA | 5 | clinical | 16 |
| Full | 5 | rna | 200 |
| Full | 5 | mirna | 7 |
| Full | 5 | methylation | 200 |
| Full | 5 | clinical | 16 |
| Clinical + Methylation | 5 | methylation | 200 |
| Clinical + Methylation | 5 | clinical | 16 |
| Clinical + miRNA | 5 | mirna | 7 |
| Clinical + miRNA | 5 | clinical | 16 |

Methylation-containing fold/configuration rows with zero selected probes: 0.
Methylation-containing fold/configuration rows with fewer than five selected probes: 0.
