# DeepSurv internal-split and feature-selection audit

## Scope and reproducibility

- Inspected code: `src/common/phase1_common.py`, `src/models/ablation_study.py`, `src/models/train_explainable_dl.py`, and `src/evaluation/cross_validate.py`.
- Outer folds: five-fold `StratifiedGroupKFold`, event-stratified and grouped by 12-character participant ID, shuffled with seed 42.
- DeepSurv internal split: `train_test_split(test_size=0.15, stratify=event, random_state=42 + outer_fold)` in the ablation pipeline.

## Verified execution order

1. Each outer fold is split into outer-training and untouched outer-test records.
2. Methylation variance filtering/imputation, clinical encoding, and univariate Cox ranking are fitted on the **entire outer-training fold**.
3. The selected outer-training matrix is passed to the DeepSurv trainer.
4. Only then is that matrix split into the neural-network fit subset (85%) and internal early-stopping validation subset (15%).
5. Neural-network scaling is fitted on the 85% fit subset only; early stopping uses the internal validation C-index.

## Finding

The survival outcomes of the later internal-validation records are included in the univariate Cox feature ranking because feature selection occurs before the internal 85/15 split. The prior statement that validation outcomes were not used by the ranking was incorrect. This is **minor leakage into the internal early-stopping validation assessment**: the validation set is not independent of supervised feature selection.

The untouched outer-test fold is excluded from variance filtering, imputation fitting, clinical encoding, univariate ranking, neural-network scaling, training, and early stopping. Therefore this issue does **not directly leak outer-test outcomes or features into model fitting**, and the reported outer-fold test C-indices remain out-of-sample estimates for the complete outer-training procedure. It can nevertheless make the internal early-stopping signal optimistic and thereby indirectly affect the epoch selected for evaluation on the outer test fold; its numerical magnitude cannot be inferred without a corrected rerun.

## Recommended correction (not executed in Phase 5)

Create the event-stratified, participant-grouped internal split before supervised feature selection; fit ranking and all inner preprocessing using only the inner-training subset; transform/select the internal-validation and outer-test subsets with those fitted objects. A corrected full rerun would be required to quantify the effect. Per the Phase 5 instruction, the model code and downstream results were not changed here.
