# SHAP v2 — Primary-tumor-only

Explainer: `shap.GradientExplainer`. Background and explained cases were random without replacement with logged fold seeds.
Mean pairwise top-20 SHAP Jaccard: 0.023627676259255204.

Feature-selection Jaccard and SHAP-ranking Jaccard are different quantities; no equivalence is inferred. Block totals and per-feature means are both reported in JSON to expose feature-count effects.

| Rank | Feature | Mean | Folds selected |
|---:|---|---:|---:|
| 1 | `rna_ENSG00000233304.7` | 0.3171467802646922 | 1 |
| 2 | `meth_cg10973622` | 0.14869974842750774 | 1 |
| 3 | `meth_cg22001496` | 0.0803323919643539 | 1 |
| 4 | `meth_cg04985582` | 0.06011568231953613 | 1 |
| 5 | `mirna_hsa-mir-7108` | 0.05534146782167578 | 1 |
| 6 | `mirna_hsa-mir-302b` | 0.03941075198342298 | 1 |
| 7 | `rna_ENSG00000250765.6` | 0.03829092100972675 | 2 |
| 8 | `rna_ENSG00000228361.2` | 0.037434845355414165 | 1 |
| 9 | `rna_ENSG00000122025.15` | 0.03737127567248217 | 1 |
| 10 | `rna_ENSG00000253105.6` | 0.037313559094074966 | 1 |
| 11 | `meth_cg23866687` | 0.03672577141049922 | 2 |
| 12 | `rna_ENSG00000072080.11` | 0.03616340450933983 | 1 |
| 13 | `rna_ENSG00000285715.1` | 0.03602330499214189 | 1 |
| 14 | `rna_ENSG00000183255.12` | 0.03601249888654153 | 1 |
| 15 | `mirna_hsa-mir-7107` | 0.03587230619554666 | 1 |
| 16 | `meth_cg07972135` | 0.031006463036382197 | 2 |
| 17 | `rna_ENSG00000229145.1` | 0.02929790699492405 | 1 |
| 18 | `meth_cg04747226` | 0.02848186778585429 | 3 |
| 19 | `mirna_hsa-mir-6846` | 0.02720725520907544 | 1 |
| 20 | `rna_ENSG00000076043.11` | 0.027159261508586323 | 1 |
