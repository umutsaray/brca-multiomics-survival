# CoxPH Convergence Diagnostics

Completed nested full-omics tuner telemetry plus clinical-only outer design audits.

## Retained-telemetry limitation

The completed tuner did not retain Newton iteration counts or a row for every inner fit. These cannot be recovered after fitting, so the report records them as unavailable rather than estimating them. A future rerun must capture `show_progress` output per fit.

## mixed: full omics nested CoxPH

| Outer fold | Selected penalizer | Warning | Iterations |
|---:|---:|---|---|
| 1 | 1.0 | True | unavailable |
| 2 | 0.5 | True | unavailable |
| 3 | 1.0 | False | unavailable |
| 4 | 0.5 | True | unavailable |
| 5 | 0.5 | True | unavailable |

Inner+outer warning counts retained by penalizer: `{'0.01': 3, '0.05': 4, '0.1': 6, '0.5': 12, '1.0': 6}`.

## mixed: clinical-only design matrix

| Fold | Rows | Columns | Rank | Rank deficient | High-correlation pairs |
|---:|---:|---:|---:|---|---:|
| 1 | 674 | 15 | 12 | True | 1 |
| 2 | 674 | 16 | 13 | True | 1 |
| 3 | 674 | 16 | 13 | True | 1 |
| 4 | 673 | 16 | 13 | True | 1 |
| 5 | 673 | 16 | 13 | True | 1 |

## primary: full omics nested CoxPH

| Outer fold | Selected penalizer | Warning | Iterations |
|---:|---:|---|---|
| 1 | 1.0 | True | unavailable |
| 2 | 1.0 | True | unavailable |
| 3 | 0.5 | True | unavailable |
| 4 | 0.5 | True | unavailable |
| 5 | 0.5 | True | unavailable |

Inner+outer warning counts retained by penalizer: `{'0.01': 2, '0.05': 8, '0.1': 5, '0.5': 13, '1.0': 6}`.

## primary: clinical-only design matrix

| Fold | Rows | Columns | Rank | Rank deficient | High-correlation pairs |
|---:|---:|---:|---:|---|---:|
| 1 | 609 | 16 | 13 | True | 1 |
| 2 | 609 | 15 | 12 | True | 1 |
| 3 | 610 | 16 | 13 | True | 1 |
| 4 | 610 | 16 | 13 | True | 1 |
| 5 | 610 | 16 | 13 | True | 1 |

## Elastic-net Cox sensitivity

Rank deficiency was confirmed in every clinical design audit. A clinical-only `CoxnetSurvivalAnalysis` sensitivity used `l1_ratio=0.5`, `alpha=0.01`, fold-local scaling, and the same outer folds. Exact results: `{'mixed': {'l1_ratio': 0.5, 'alpha': 0.01, 'fold_cindices': [0.6723735408560312, 0.7526835551738944, 0.7199275362318841, 0.7986735130509199, 0.7497123130034522], 'mean': 0.7386740916632364, 'sample_sd': 0.04654362024743143}, 'primary': {'l1_ratio': 0.5, 'alpha': 0.01, 'fold_cindices': [0.8230717185385656, 0.6186813186813187, 0.8103693181818182, 0.7371883656509696, 0.80375], 'mean': 0.7586121442105345, 'sample_sd': 0.0850160014737996}}`.

## Rare source categories

Exact n<10 source-level counts are stored in `coxph_convergence_diagnostics.json`.