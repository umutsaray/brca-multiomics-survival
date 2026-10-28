# Phase 16 Task 0(a) — HER2 source cross-tabulation

## Pre-specified reporting rules

- Rows: `lab_proc_her2_neu_immunohistochemistry_receptor_status`; columns: `lab_procedure_her2_neu_in_situ_hybrid_outcome_type`.
- `blank/NA` is explicit. No values were recoded and no model was rerun.
- Distinct normalized nonblank categories are labelled `genuine category disagreement`; identical normalized categories are `same normalized category / wording only`.
- The Phase 15 conflict rule used three sources, so the 163 participants are decomposed by the full IHC/ISH/nature2012 tuple below.

## IHC × ISH cross-tabulation (all 762 participants)

| IHC \ ISH | `Equivocal` | `Indeterminate` | `Negative` | `Positive` | `blank/NA` | Row total |
|---|---|---|---|---|---|---|
| `Equivocal` | 3 | 0 | 106 | 16 | 14 | 139 |
| `Indeterminate` | 0 | 0 | 6 | 2 | 4 | 12 |
| `Negative` | 0 | 3 | 49 | 2 | 331 | 385 |
| `Positive` | 1 | 1 | 8 | 10 | 71 | 91 |
| `blank/NA` | 0 | 0 | 45 | 6 | 84 | 135 |
| Column total | 4 | 4 | 214 | 36 | 504 | 762 |

## Phase 15 HER2 Conflict participants by full value tuple

| IHC | ISH | HER2 Final nature2012 | n | classification |
|---|---|---|---|---|
| `Equivocal` | `Negative` | `Negative` | 70 | genuine category disagreement |
| `Equivocal` | `Negative` | `blank/NA` | 36 | genuine category disagreement |
| `Equivocal` | `Positive` | `Positive` | 10 | genuine category disagreement |
| `Positive` | `blank/NA` | `Negative` | 10 | genuine category disagreement |
| `Equivocal` | `Positive` | `blank/NA` | 6 | genuine category disagreement |
| `Positive` | `Negative` | `Negative` | 6 | genuine category disagreement |
| `Indeterminate` | `Negative` | `Negative` | 5 | genuine category disagreement |
| `Equivocal` | `blank/NA` | `Negative` | 3 | genuine category disagreement |
| `Equivocal` | `blank/NA` | `Positive` | 3 | genuine category disagreement |
| `Negative` | `Indeterminate` | `blank/NA` | 3 | genuine category disagreement |
| `Indeterminate` | `Positive` | `Positive` | 2 | genuine category disagreement |
| `Equivocal` | `Equivocal` | `Negative` | 1 | genuine category disagreement |
| `Indeterminate` | `Negative` | `blank/NA` | 1 | genuine category disagreement |
| `Negative` | `Positive` | `Negative` | 1 | genuine category disagreement |
| `Negative` | `Positive` | `blank/NA` | 1 | genuine category disagreement |
| `Positive` | `Equivocal` | `blank/NA` | 1 | genuine category disagreement |
| `Positive` | `Indeterminate` | `blank/NA` | 1 | genuine category disagreement |
| `Positive` | `Negative` | `Positive` | 1 | genuine category disagreement |
| `Positive` | `Negative` | `blank/NA` | 1 | genuine category disagreement |
| `Positive` | `Positive` | `Negative` | 1 | genuine category disagreement |

Tuple count sum: `163`; Phase 15 Conflict count: `163`.

## Reproducibility

Command: `src\evaluation\phase16_task0_missing_phase15_items.py --repo-dir .`  
Versions: `{"python": "3.14.7", "pandas": "2.3.3"}`  
Seeds: `null`.  
Runtime seconds: `0.0567548000253737`.
