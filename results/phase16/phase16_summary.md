# Phase 16 summary

Values are copied from Phase 16 JSON outputs without display rounding.

## Task 0
- Status: completed; HER2 conflicts=163; full tuple count sum=163.

## Task 3
- Status: completed; runtime_seconds=1110.4178037999664.
- Late fusion — equal Cox risks, coxph: Clinical-only=0.7598208476599109; new=0.7092179551332813; Clinical-only-new=0.05060289252662953; ordering changed=no.
- Late fusion — stacked Cox risks, coxph: Clinical-only=0.7598208476599109; new=0.7028794600559853; Clinical-only-new=0.05694138760392553; ordering changed=no.
- PCA early fusion k=5, coxph: Clinical-only=0.7598208476599109; new=0.7128704242431907; Clinical-only-new=0.04695042341672018; ordering changed=no.
- PCA early fusion k=5, rsf: Clinical-only=0.7572530299391804; new=0.6065568593798887; Clinical-only-new=0.15069617055929174; ordering changed=no.
- PCA early fusion k=10, coxph: Clinical-only=0.7598208476599109; new=0.6916574720625898; Clinical-only-new=0.06816337559732111; ordering changed=no.
- PCA early fusion k=10, rsf: Clinical-only=0.7572530299391804; new=0.5728098878644573; Clinical-only-new=0.18444314207472312; ordering changed=no.
- PCA early fusion k=20, coxph: Clinical-only=0.7598208476599109; new=0.6900609062620214; Clinical-only-new=0.06975994139788944; ordering changed=no.
- PCA early fusion k=20, rsf: Clinical-only=0.7572530299391804; new=0.587749799576811; Clinical-only-new=0.16950323036236947; ordering changed=no.

## Task 4
- Status: partially completed; runtime_seconds=1486.2478894998785; stopped={'Variant E': 'Stopped before any outer-fold result was produced; projected runtime exceeded the approximately 4-hour notification/feasibility threshold after the first 6000-feature elastic-net inner path did not finish in approximately 12 minutes. No substitute analysis was used.', 'Variant V DeepSurv': 'Skipped for runtime as pre-specified.'}.
- Variant V — variance top200/block, coxph: Clinical-only=0.7598208476599109; new=0.699078990335383; Clinical-only-new=0.06074185732452786; ordering changed=no.
- Variant V — variance top200/block, rsf: Clinical-only=0.7572530299391804; new=0.677297503600965; Clinical-only-new=0.07995552633821545; ordering changed=no.

## Task 5
- Status: completed; runtime_seconds=22.257798000006005; stored C-index reproduced exactly=True.
- Fold 1: within-fold top-20 Jaccard=0.47160535049257357; rank Spearman=0.6698647986969513.
- Fold 2: within-fold top-20 Jaccard=0.2524134872463081; rank Spearman=0.5300584351201401.
- Fold 3: within-fold top-20 Jaccard=0.3449137335059965; rank Spearman=0.47635080805668334.
- Fold 4: within-fold top-20 Jaccard=0.2493438903922541; rank Spearman=0.36823032946312506.
- Fold 5: within-fold top-20 Jaccard=0.29447146673512936; rank Spearman=0.4720863005235261.
- Selected in all five folds=15; cross-fold common-feature rank Spearman=0.6649999999999998.
