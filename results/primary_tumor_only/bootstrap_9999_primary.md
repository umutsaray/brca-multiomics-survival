# 9,999-iteration paired participant bootstrap

- Iterations requested: 9999
- Random seed: 42
- Resampling unit: 12-character TCGA participant ID
- Comparison direction: Full multi-omics minus Clinical-only
- Input risks: unchanged stored out-of-fold predictions

## mixed

| Model | Δ C-index | 95% CI | Two-sided bootstrap p-value | Valid iterations |
|---|---:|---|---:|---:|
| coxph | -0.0681665036409086 | [-0.13728548482841296, -0.0019052811886452984] | 0.0448 | 9999 |
| rsf | -0.06089555483099662 | [-0.1308877898320284, 0.010998422960275011] | 0.095 | 9999 |
| deepsurv | -0.01613955004890777 | [-0.12002255936521569, 0.0835761367061322] | 0.7692 | 9999 |

## primary

| Model | Δ C-index | 95% CI | Two-sided bootstrap p-value | Valid iterations |
|---|---:|---|---:|---:|
| coxph | -0.12286610685614097 | [-0.1884310470168456, -0.05496531706757675] | 0.0002 | 9999 |
| rsf | -0.14727323059887432 | [-0.21664871667470406, -0.07350564789056625] | 0.0004 | 9999 |
| deepsurv | -0.214850358340254 | [-0.27757620412555767, -0.15179862417960882] | 0.0002 | 9999 |
