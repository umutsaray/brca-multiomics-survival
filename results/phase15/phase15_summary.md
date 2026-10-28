# Phase 15 summary

All values below are copied or deterministically subtracted from the Phase 15 JSON outputs; no display rounding was applied.

## Task 1 — ER/PR/HER2 missingness

- Status: completed.
- Matched corrected primary-tumor-only participants: 762 / 762.
- ER: all participants={'Positive': 560, 'Negative': 162, 'Not available': 40}; events={'Positive': 61, 'Negative': 31, 'Not available': 7}; conflicts=0.
- PR: all participants={'Positive': 490, 'Negative': 229, 'Not available': 43}; events={'Positive': 54, 'Negative': 37, 'Not available': 8}; conflicts=0.
- HER2: all participants={'Positive': 80, 'Negative': 479, 'Not available': 203}; events={'Positive': 10, 'Negative': 65, 'Not available': 24}; conflicts=163.
- Complete-case cohort: n=530; events=72.
- PAM50 used as model input: False.

## Task 2A — Missing-indicator analysis (all participants)

- Status: completed; n=762; events=99; runtime_seconds=638.6391574000008.
- Table 7 reproduction, Clinical-only, coxph: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- Table 7 reproduction, Clinical-only, rsf: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- Table 7 reproduction, Clinical-only, deepsurv: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- Table 7 reproduction, Full, coxph: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- Table 7 reproduction, Full, rsf: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- Table 7 reproduction, Full, deepsurv: exact_match=True; differences=[0.0, 0.0, 0.0, 0.0, 0.0].
- 2A, coxph: Clinical-only=0.7598208476599109; Full=0.6524726851411089; Clinical-only - Full=0.10734816251880197; ordering changed=no (this is the unchanged main-analysis ordering).
- 2A, rsf: Clinical-only=0.7572530299391804; Full=0.6220342521048454; Clinical-only - Full=0.135218777834335; ordering changed=no (this is the unchanged main-analysis ordering).
- 2A, deepsurv: Clinical-only=0.7443252097099942; Full=0.5262759610444295; Clinical-only - Full=0.21804924866556463; ordering changed=no (this is the unchanged main-analysis ordering).

## Task 2B — Complete-case analysis

- Status: completed; n=530; events=72; runtime_seconds=6310.380588200176.
- Outer-test-fold event counts: [14, 15, 15, 14, 14]; folds with fewer than 5 events: [].
- 2B (relative to 2A main ordering), coxph: Clinical-only=0.7368866149723537; Full=0.6055094029868497; Clinical-only - Full=0.131377211985504; ordering changed=no.
- 2B (relative to 2A main ordering), rsf: Clinical-only=0.7077419855122589; Full=0.6080089585487123; Clinical-only - Full=0.09973302696354658; ordering changed=no.
- 2B (relative to 2A main ordering), deepsurv: Clinical-only=0.6470473735705957; Full=0.5511203940914301; Clinical-only - Full=0.09592697947916562; ordering changed=no.

## Files

- `task1_receptor_missingness.md` / `.json`
- `task2a_missing_indicator.md` / `.json`
- `task2b_complete_case.md` / `.json`
