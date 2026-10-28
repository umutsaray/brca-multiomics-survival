# Phase 17 Task 6 — Reduced CoxPH-only learning curve

Status: **stopped**.

RSF and DeepSurv were skipped as a pre-fixed runtime decision. Fractions 50% and 75% use three repetitions; 100% is the stored Table 7 CoxPH result.

The reduced design did not complete. Available reduced-fraction values are checkpointed partial results and are not interpreted as the completed learning-curve estimand.

Stop reason: Stopped because projected remaining runtime 30490.219383898893 seconds exceeded the pre-fixed 25200.0-second threshold.
Projected remaining seconds at stop: `30490.219383898893`.

## C-index and training-size summary

| Fraction | Configuration | Completed folds | Mean C-index | Sample SD | Mean training participants/fold | Mean training events/fold | Mean selected features/block | Fold values |
|---:|---|---:|---:|---:|---:|---:|---|---|
| 0.5 | Clinical-only | 1 | 0.7822507893549843 | None | 304.0 | 39.0 | `{"rna": 0.0, "mirna": 0.0, "methylation": 0.0, "clinical": 15.0}` | `[0.7822507893549843]` |
| 0.5 | Full | 1 | 0.6783942264321153 | None | 304.0 | 39.0 | `{"rna": 160.66666666666666, "mirna": 6.666666666666667, "methylation": 25.0, "clinical": 15.0}` | `[0.6783942264321153]` |
| 0.75 | Clinical-only | 0 | None | None | None | None | `null` | `[]` |
| 0.75 | Full | 0 | None | None | None | None | `null` | `[]` |
| 1.0 | Clinical-only | 5 | 0.7598208476599109 | 0.07678272079706581 | 609.6 | 79.2 | `{"rna": 0.0, "mirna": 0.0, "methylation": 0.0, "clinical": 15.8}` | `[0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375]` |
| 1.0 | Full | 5 | 0.6524726851411089 | 0.06777572101415535 | 609.6 | 79.2 | `{"rna": 181.4, "mirna": 8.6, "methylation": 3.4, "clinical": 15.8}` | `[0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667]` |

## Clinical-only minus Full bootstrap

| Fraction | Point difference | 95% CI |
|---:|---:|---|
| 0.5 | None | `null` |
| 0.75 | None | `null` |
| 1.0 | 0.10734816251880197 | `[0.03700289487664033, 0.17268175928129118]` |

## Observed trends

- Clinical-only: `None`; means at 50%, 75%, 100% = `[0.7822507893549843, null, 0.7598208476599109]`.
- Full: `None`; means at 50%, 75%, 100% = `[0.6783942264321153, null, 0.6524726851411089]`.

Per-fold training participant/event counts, repetition-level C-indices, selected feature counts per block, raw risks, and all bootstrap draws are retained without rounding in the JSON/checkpoints.

Actual elapsed time (seconds): `3525.293533799937`.
