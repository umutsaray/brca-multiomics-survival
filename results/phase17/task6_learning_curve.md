# Phase 17 Task 6 — Reduced CoxPH-only learning curve

Status: **stopped**.

The 50% values contain one completed outer fold and one repetition only, not the pre-specified five-fold repetition-averaged estimand; the 75% analysis was not started. Neither is interpreted as a completed learning-curve result.

RSF and DeepSurv were skipped as a pre-fixed runtime decision. CoxPH used the Table 7 fixed ridge penalizer 0.1. Fractions 50% and 75% were specified for three repetitions; 100% is the stored Table 7 result.

Stop reason: Stopped because projected remaining runtime 20513.152537505142 seconds exceeded the pre-fixed 10800-second threshold.
Projected remaining seconds at stop: `20513.152537505142`.
Completed checkpoint count: `1`.

## C-index summary

| Fraction | Configuration | Completed folds | Mean | Sample SD | Fold values |
|---:|---|---:|---:|---:|---|
| 0.5 | Clinical-only | 1 | 0.770297699594046 | None | `[0.770297699594046]` |
| 0.5 | Full | 1 | 0.7253044654939107 | None | `[0.7253044654939107]` |
| 0.75 | Clinical-only | 0 | None | None | `[]` |
| 0.75 | Full | 0 | None | None | `[]` |
| 1.0 | Clinical-only | 5 | 0.7598208476599109 | 0.07678272079706581 | `[0.8271312584573748, 0.6362637362637362, 0.7958096590909091, 0.7361495844875346, 0.80375]` |
| 1.0 | Full | 5 | 0.6524726851411089 | 0.06777572101415535 | `[0.7401894451962111, 0.6615384615384615, 0.6200284090909091, 0.6814404432132964, 0.5591666666666667]` |

## Completed reduced-fit checkpoint

| Fraction | Fold | Repetition | Configuration | Training participants | Training events | Selected RNA | Selected miRNA | Selected methylation | Selected clinical | C-index |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 1 | 1 | Clinical-only | 304 | 39 | 0 | 0 | 0 | 16 | 0.770297699594046 |
| 0.5 | 1 | 1 | Full | 304 | 39 | 200 | 9 | 75 | 16 | 0.7253044654939107 |

## Clinical-only minus Full fold-stratified paired bootstrap

| Fraction | Point difference | 95% CI | Status |
|---:|---:|---|---|
| 0.5 | None | `null` | unavailable: incomplete design |
| 0.75 | None | `null` | unavailable: incomplete design |
| 1.0 | 0.10734816251880197 | `[0.03700289487664033, 0.17268175928129118]` | available |

## Observed trend

No 50%–75%–100% trend classification is reported because the 50% design is incomplete and the 75% design was not started. No extrapolation was performed.

Actual elapsed time (seconds): `838.0257824000437`.

Report finalization: `{"command": "C:\\moenv\\Scripts\\python.exe src/evaluation/phase17_task6_finalize.py --json results/phase17/task6_learning_curve.json --markdown results/phase17/task6_learning_curve.md --checkpoint-dir results/phase17/task6_checkpoints", "model_fits_run": 0, "bootstrap_resamples_run": 0, "derived_fields": ["correct checkpoint-file count", "stored-summary Clinical-only minus Full point difference"], "statistics_replaced": false}`.
