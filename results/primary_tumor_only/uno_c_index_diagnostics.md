# Uno C-index diagnostics

- Script: `src/evaluation/uno_c_index_diagnostics.py`
- Inputs: retained outer-fold OOF predictions from `ablation_oof_extended.parquet` and `ablation_oof_primary.parquet`.
- Fold definition: the stored five participant-grouped outer folds (global split seed 42).
- The retained analysis called `concordance_index_ipcw` with `tau=None`, so no fixed truncation horizon was imposed. `effective_max_event_time_days` is the latest held-fold event time that can contribute; maximum observed follow-up is reported separately.
- IPCW values are estimated from the other four folds; zero entries correspond to censored observations. The concordance implementation applies squared IPCW event weights.

| cohort | configuration | fold | n_test | effective_max_event_time_days | max_observed_time_days | at_risk_at_effective_horizon | events_by_effective_horizon | event_ipcw_min | event_ipcw_median | event_ipcw_max | concordance_weight_max_ipcw_squared | n_nonzero_event_weights |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mixed | Clinical-only | 1 | 168 | 6593.0 | 6593.0 | 1 | 26 | 1.0634041766532523 | 4.040005552858319 | 68.41985879361502 | 4681.277077338219 | 26 |
| mixed | Clinical-only | 2 | 168 | 2712.0 | 7777.0 | 17 | 26 | 1.0014858841010401 | 1.9553257051807384 | 6.075954130514331 | 36.91721859611416 | 26 |
| mixed | Clinical-only | 3 | 168 | 3959.0 | 4894.0 | 4 | 26 | 1.1348231427296325 | 2.4578517591906723 | 19.188522802335733 | 368.19940733575834 | 26 |
| mixed | Clinical-only | 4 | 169 | 3669.0 | 8008.0 | 7 | 27 | 1.04031007751938 | 1.9751041275111756 | 16.99106994263566 | 288.696457795537 | 27 |
| mixed | Clinical-only | 5 | 169 | 3262.0 | 8605.0 | 6 | 27 | 1.0484629621112873 | 2.518570011522922 | 10.029883609069955 | 100.59856521149015 | 27 |
| mixed | Full | 1 | 168 | 6593.0 | 6593.0 | 1 | 26 | 1.0634041766532523 | 4.040005552858319 | 68.41985879361502 | 4681.277077338219 | 26 |
| mixed | Full | 2 | 168 | 2712.0 | 7777.0 | 17 | 26 | 1.0014858841010401 | 1.9553257051807384 | 6.075954130514331 | 36.91721859611416 | 26 |
| mixed | Full | 3 | 168 | 3959.0 | 4894.0 | 4 | 26 | 1.1348231427296325 | 2.4578517591906723 | 19.188522802335733 | 368.19940733575834 | 26 |
| mixed | Full | 4 | 169 | 3669.0 | 8008.0 | 7 | 27 | 1.04031007751938 | 1.9751041275111756 | 16.99106994263566 | 288.696457795537 | 27 |
| mixed | Full | 5 | 169 | 3262.0 | 8605.0 | 6 | 27 | 1.0484629621112873 | 2.518570011522922 | 10.029883609069955 | 100.59856521149015 | 27 |
| primary | Clinical-only | 1 | 153 | 6593.0 | 7777.0 | 2 | 20 | 1.032258064516129 | 1.734747193917681 | 100.6695433324238 | 10134.356954758752 | 20 |
| primary | Clinical-only | 2 | 153 | 3669.0 | 5739.0 | 4 | 20 | 1.0611652235258158 | 5.089874147597087 | 20.328121610833378 | 413.232528224831 | 20 |
| primary | Clinical-only | 3 | 152 | 3262.0 | 8605.0 | 6 | 20 | 1.0339618080734831 | 2.039316326818023 | 10.583694668845325 | 112.01459284334496 | 20 |
| primary | Clinical-only | 4 | 152 | 2192.0 | 5156.0 | 28 | 19 | 1.0628657801843584 | 2.2092236697554197 | 4.491375991130539 | 20.172458293703833 | 19 |
| primary | Clinical-only | 5 | 152 | 3959.0 | 8008.0 | 5 | 20 | 1.0032894736842104 | 2.817462357374378 | 22.37789990070304 | 500.77040396588507 | 20 |
| primary | Full | 1 | 153 | 6593.0 | 7777.0 | 2 | 20 | 1.032258064516129 | 1.734747193917681 | 100.6695433324238 | 10134.356954758752 | 20 |
| primary | Full | 2 | 153 | 3669.0 | 5739.0 | 4 | 20 | 1.0611652235258158 | 5.089874147597087 | 20.328121610833378 | 413.232528224831 | 20 |
| primary | Full | 3 | 152 | 3262.0 | 8605.0 | 6 | 20 | 1.0339618080734831 | 2.039316326818023 | 10.583694668845325 | 112.01459284334496 | 20 |
| primary | Full | 4 | 152 | 2192.0 | 5156.0 | 28 | 19 | 1.0628657801843584 | 2.2092236697554197 | 4.491375991130539 | 20.172458293703833 | 19 |
| primary | Full | 5 | 152 | 3959.0 | 8008.0 | 5 | 20 | 1.0032894736842104 | 2.817462357374378 | 22.37789990070304 | 500.77040396588507 | 20 |

## Observed weight ranges

| cohort | configuration | max_event_ipcw | max_squared_weight |
|---|---|---|---|
| mixed | Clinical-only | 68.41985879361502 | 4681.277077338219 |
| mixed | Full | 68.41985879361502 | 4681.277077338219 |
| primary | Clinical-only | 100.6695433324238 | 10134.356954758752 |
| primary | Full | 100.6695433324238 | 10134.356954758752 |

Large weights are reported descriptively only. Association with the Harrell–Uno discrepancy is not asserted without a dedicated perturbation analysis.
