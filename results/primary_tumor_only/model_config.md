# Model and Validation Configuration

This document records the settings actually used by the Phase 1/2 code. Global CLI default seed is 42. Outer splitting is five-fold `StratifiedGroupKFold(shuffle=True, random_state=42)`, stratified on event and grouped by the 12-character TCGA participant ID. In the primary-only cohort there is exactly one record per participant.

## CoxPH

- Implementation: `lifelines.CoxPHFitter`.
- Fixed-penalty ablation/sensitivity runs: ridge penalizer 0.1, `l1_ratio=0`, maximum Newton steps 50.
- Nested tuning candidates: 0.01, 0.05, 0.1, 0.5, and 1.0. Inner CV is three-fold `StratifiedGroupKFold`, grouped by participant, shuffled with seed `42 + outer_fold`.
- Reported risk is `predict_partial_hazard`, not log-risk.
- Continuous predictors are standardized using the corresponding fitting data only.
- Important limitation: in the current nested tuner, methylation variance filtering/imputation, univariate feature selection, clinical encoding, and the outer-fold feature list are fitted once on the outer training set. They are not re-fitted inside each inner split. Consequently, the inner penalizer-selection score is not a fully nested estimate for preprocessing/feature selection; the untouched outer test score remains the evaluation target. This limitation must be disclosed.

## Random Survival Forest

- Implementation: `sksurv.ensemble.RandomSurvivalForest`.
- `n_estimators=500`, `min_samples_split=10`, `min_samples_leaf=5`, `max_features='sqrt'` (library default), `max_depth=None`, bootstrap enabled, `n_jobs=-1`.
- `random_state=42 + outer_fold` for ablation; feature-count sensitivity uses `42 + 100*outer_fold + N`.
- Reported risk is the estimator's `predict` risk score.

## DeepSurv-style neural network

- Separate two-layer encoders: RNA 128→32, miRNA 64→16, methylation 128→32, clinical 32→8. Each encoder layer is Linear–BatchNorm–ReLU–Dropout.
- Concatenated head: Linear(total, 64)–ReLU–Dropout–Linear(64, 16)–ReLU–Linear(16, 1).
- Dropout 0.3. Optimizer AdamW, learning rate 0.001, weight decay 0.0001.
- Maximum 300 epochs; early-stopping patience 30; improvement tolerance 0.000001.
- Training is full-batch: batch size equals the fit subset size. The risk set therefore spans the complete fit subset at every optimization step.
- Loss is the mean negative Cox partial log-likelihood. Event times are sorted descending; `logcumsumexp` forms risk sets. Ties use the code's Breslow-style shared full tied risk set.
- Internal validation fraction is 0.15 and is event-stratified via `train_test_split`. It is not explicitly participant-grouped. This is equivalent to participant splitting only in the primary-only cohort; in the mixed cohort, duplicate records from a participant can cross this internal split.
- Univariate selection and outer-fold imputation are performed before this internal validation split. Therefore, the survival outcomes of records later assigned to internal validation are used by the univariate Cox ranking. This creates minor leakage into the internal early-stopping validation assessment, although the untouched outer-test fold remains excluded from preprocessing, selection, training, and early stopping. Scaling is fitted only on the neural-network fit subset and then applied to validation/test. See `deepsurv_internal_split_audit.md`.
- Seeds set Python `random`, NumPy, and PyTorch; deterministic algorithms are requested with warning-only fallback. Ablation seed is `42 + outer_fold`; sensitivity seed is `42 + 100*outer_fold + N`.
- Reported risk is the network's raw scalar output (log-risk scale).

## Fold-local versus inner-fold preprocessing

For every outer fold, methylation variance selection, KNN imputation, clinical imputation/one-hot encoding, univariate Cox selection, and scaling are learned without the outer test fold. In the nested Cox penalizer tuner, these steps are not repeated inside the inner folds, as stated above. There is no inner hyperparameter CV for RSF or DeepSurv. Clinical AJCC stages are collapsed to Stage I/II/III/IV after printing the categories actually present in the source file.

## Feature-selection procedures and the former Table 3b/Table 4 mismatch

`select_features_v2.py` is a feature-count sensitivity analysis: it takes raw-p-value-ranked top N features per omics block without an FDR cutoff and changes the model seed by N. `ablation_study.py` uses BH-FDR < 0.05 and then caps each block at 200 features, with a different model seed schedule. Both use the same outer participant folds, but they do not estimate the same feature set/model realization. Their N=200/full values therefore must not be presented as duplicate estimates. The FDR ablation is the canonical full-model estimate for comparisons with clinical-only; the raw-p top-N result is a sensitivity analysis.
