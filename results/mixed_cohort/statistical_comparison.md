# Statistical Comparison of Survival Models

The Wilcoxon signed-rank tests use only five paired CV folds (`n=5`). Therefore statistical power is low and non-significance must not be interpreted as equivalence.
Bootstrap comparisons use the same 1,000 patient resample indices for all three models.

| Comparison | Wilcoxon W | Wilcoxon p | Observed test ΔC | Bootstrap mean ΔC | Paired-bootstrap 95% CI | Empirical bootstrap p | Conclusion |
|---|---:|---:|---:|---:|---:|---:|---|
| DeepSurv vs CoxPH | 3 | 0.3125 | +0.039216 | +0.039348 | [-0.053161, +0.166321] | 0.521479 | Anlamlı fark bulunamadı; paired bootstrap %95 CI sıfırı içeriyor. |
| DeepSurv vs Random Survival Forest | 1 | 0.125 | +0.030897 | +0.031312 | [-0.064761, +0.125915] | 0.495504 | Anlamlı fark bulunamadı; paired bootstrap %95 CI sıfırı içeriyor. |

## Interpretation

- **DeepSurv vs CoxPH:** Anlamlı fark bulunamadı; paired bootstrap %95 CI sıfırı içeriyor.
- **DeepSurv vs Random Survival Forest:** Anlamlı fark bulunamadı; paired bootstrap %95 CI sıfırı içeriyor.
