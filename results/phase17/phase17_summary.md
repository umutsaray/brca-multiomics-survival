# Phase 17 summary

## Task 7 — Uno versus Harrell sensitivity

Status: **completed**. Stored OOF predictions only; no model refitting.

Runtime seconds: `8039.219540800201`. Bootstrap: `9999` requested iterations, seed `42`.

Each line below states whether the Clinical-only versus Full ordering changed relative to the existing-tau Uno sign under the same cohort/model/weight construction.

| Cohort | Model | Tau | Weight construction | Clinical-only − Full | 95% CI | Valid iterations | Ordering changed? |
|---|---|---|---|---:|---|---:|---|
| mixed | coxph | existing_none | other_four_oof | -0.019353720696004207 | `[-0.10509036106493296, 0.0791708482505417]` | 9999 | no |
| mixed | coxph | existing_none | outer_training_outcomes | -0.019353720696004207 | `[-0.10509036106493296, 0.0791708482505417]` | 9999 | no |
| mixed | coxph | 3_year | other_four_oof | 0.008741796960769044 | `[-0.08014464719511778, 0.10670251350061215]` | 9972 | yes |
| mixed | coxph | 3_year | outer_training_outcomes | 0.008741796960769044 | `[-0.08014464719511778, 0.10670251350061215]` | 9972 | yes |
| mixed | coxph | 5_year | other_four_oof | 0.06563693517110736 | `[-0.020146286451139857, 0.15395506893465352]` | 9998 | yes |
| mixed | coxph | 5_year | outer_training_outcomes | 0.06563693517110736 | `[-0.020146286451139857, 0.15395506893465352]` | 9998 | yes |
| mixed | coxph | 7_year | other_four_oof | 0.04061283482310864 | `[-0.036221445121145536, 0.1229057158675549]` | 9999 | yes |
| mixed | coxph | 7_year | outer_training_outcomes | 0.04061283482310864 | `[-0.036221445121145536, 0.1229057158675549]` | 9999 | yes |
| mixed | rsf | existing_none | other_four_oof | -0.026936912174609473 | `[-0.12083535227712465, 0.08069741329545473]` | 9999 | no |
| mixed | rsf | existing_none | outer_training_outcomes | -0.026936912174609473 | `[-0.12083535227712465, 0.08069741329545473]` | 9999 | no |
| mixed | rsf | 3_year | other_four_oof | 0.04759695884328745 | `[-0.06378237886993106, 0.1799927419235375]` | 9972 | yes |
| mixed | rsf | 3_year | outer_training_outcomes | 0.04759695884328745 | `[-0.06378237886993106, 0.1799927419235375]` | 9972 | yes |
| mixed | rsf | 5_year | other_four_oof | 0.044171620722510775 | `[-0.04721743397636417, 0.13978343621502895]` | 9998 | yes |
| mixed | rsf | 5_year | outer_training_outcomes | 0.044171620722510775 | `[-0.04721743397636417, 0.13978343621502895]` | 9998 | yes |
| mixed | rsf | 7_year | other_four_oof | -0.021467217629513757 | `[-0.09751692674385629, 0.062608988215795]` | 9999 | no |
| mixed | rsf | 7_year | outer_training_outcomes | -0.021467217629513757 | `[-0.09751692674385629, 0.062608988215795]` | 9999 | no |
| mixed | deepsurv | existing_none | other_four_oof | -0.04978505621603002 | `[-0.1427022153784877, 0.07875228941058408]` | 9999 | no |
| mixed | deepsurv | existing_none | outer_training_outcomes | -0.04978505621603002 | `[-0.1427022153784877, 0.07875228941058408]` | 9999 | no |
| mixed | deepsurv | 3_year | other_four_oof | 0.09698403969765998 | `[-0.02219020377964897, 0.23863689545199274]` | 9972 | yes |
| mixed | deepsurv | 3_year | outer_training_outcomes | 0.09698403969765998 | `[-0.02219020377964897, 0.23863689545199274]` | 9972 | yes |
| mixed | deepsurv | 5_year | other_four_oof | 0.11787429071024744 | `[0.010503247831023691, 0.2258819115691073]` | 9998 | yes |
| mixed | deepsurv | 5_year | outer_training_outcomes | 0.11787429071024744 | `[0.010503247831023691, 0.2258819115691073]` | 9998 | yes |
| mixed | deepsurv | 7_year | other_four_oof | 0.06374413032554486 | `[-0.02703193370747662, 0.1568785384278921]` | 9999 | yes |
| mixed | deepsurv | 7_year | outer_training_outcomes | 0.06374413032554486 | `[-0.02703193370747662, 0.1568785384278921]` | 9999 | yes |
| primary | coxph | existing_none | other_four_oof | 0.016605331610658126 | `[-0.06592665131201647, 0.11701822670538543]` | 9999 | no |
| primary | coxph | existing_none | outer_training_outcomes | 0.016605331610658126 | `[-0.06592665131201647, 0.11701822670538543]` | 9999 | no |
| primary | coxph | 3_year | other_four_oof | 0.10661252477381548 | `[0.011698764265270434, 0.19547502439571712]` | 9839 | no |
| primary | coxph | 3_year | outer_training_outcomes | 0.10661252477381548 | `[0.011698764265270434, 0.19547502439571712]` | 9839 | no |
| primary | coxph | 5_year | other_four_oof | 0.10498879920702286 | `[0.01308690835494376, 0.1929483348617457]` | 9997 | no |
| primary | coxph | 5_year | outer_training_outcomes | 0.10498879920702286 | `[0.01308690835494376, 0.1929483348617457]` | 9997 | no |
| primary | coxph | 7_year | other_four_oof | 0.0913834692222349 | `[0.004200550730348545, 0.17966957719642002]` | 9999 | no |
| primary | coxph | 7_year | outer_training_outcomes | 0.0913834692222349 | `[0.004200550730348545, 0.17966957719642002]` | 9999 | no |
| primary | rsf | existing_none | other_four_oof | 0.01924583849949375 | `[-0.06743312306075193, 0.127483906523009]` | 9999 | no |
| primary | rsf | existing_none | outer_training_outcomes | 0.01924583849949375 | `[-0.06743312306075193, 0.127483906523009]` | 9999 | no |
| primary | rsf | 3_year | other_four_oof | 0.11013657801320427 | `[-0.0031884624808052, 0.21412844018882207]` | 9839 | no |
| primary | rsf | 3_year | outer_training_outcomes | 0.11013657801320427 | `[-0.0031884624808052, 0.21412844018882207]` | 9839 | no |
| primary | rsf | 5_year | other_four_oof | 0.11470267849088375 | `[0.017453026506140053, 0.20688458163527607]` | 9997 | no |
| primary | rsf | 5_year | outer_training_outcomes | 0.11470267849088375 | `[0.017453026506140053, 0.20688458163527607]` | 9997 | no |
| primary | rsf | 7_year | other_four_oof | 0.07902463886636446 | `[-0.002827072096000236, 0.16428935127754404]` | 9999 | no |
| primary | rsf | 7_year | outer_training_outcomes | 0.07902463886636446 | `[-0.002827072096000236, 0.16428935127754404]` | 9999 | no |
| primary | deepsurv | existing_none | other_four_oof | 0.07519112895366087 | `[0.000420027336828892, 0.18665923484060246]` | 9999 | no |
| primary | deepsurv | existing_none | outer_training_outcomes | 0.07519112895366087 | `[0.000420027336828892, 0.18665923484060246]` | 9999 | no |
| primary | deepsurv | 3_year | other_four_oof | 0.2502031567485644 | `[0.11406675470169332, 0.3799172360990463]` | 9839 | no |
| primary | deepsurv | 3_year | outer_training_outcomes | 0.2502031567485644 | `[0.11406675470169332, 0.3799172360990463]` | 9839 | no |
| primary | deepsurv | 5_year | other_four_oof | 0.12953916098895152 | `[0.017480171538106014, 0.25496504360720923]` | 9997 | no |
| primary | deepsurv | 5_year | outer_training_outcomes | 0.12953916098895152 | `[0.017480171538106014, 0.25496504360720923]` | 9997 | no |
| primary | deepsurv | 7_year | other_four_oof | 0.11898706088532318 | `[0.03344335482553859, 0.20938108794278035]` | 9999 | no |
| primary | deepsurv | 7_year | outer_training_outcomes | 0.11898706088532318 | `[0.03344335482553859, 0.20938108794278035]` | 9999 | no |

The two censoring-outcome constructions used identical outcome sets in every fold: mixed=`True`, primary=`True`.

## Task 6 — Reduced CoxPH-only learning curve

Status: **stopped**.

Stop reason: Stopped because projected remaining runtime 20513.152537505142 seconds exceeded the pre-fixed 10800-second threshold.
Actual elapsed seconds: `838.0257824000437`. Projected remaining seconds: `20513.152537505142`. Completed checkpoints: `1`.

RSF and DeepSurv were skipped by the pre-fixed runtime decision. CoxPH used fixed ridge penalizer 0.1.

The 50% values contain one completed outer fold and one repetition only, not the pre-specified five-fold repetition-averaged estimand; the 75% analysis was not started. Neither is interpreted as a completed learning-curve result.

| Fraction | Configuration | Mean | Sample SD | Completed folds |
|---:|---|---:|---:|---:|
| 0.5 | Clinical-only | 0.770297699594046 | None | 1 |
| 0.5 | Full | 0.7253044654939107 | None | 1 |
| 0.75 | Clinical-only | None | None | 0 |
| 0.75 | Full | None | None | 0 |
| 1.0 | Clinical-only | 0.7598208476599109 | 0.07678272079706581 | 5 |
| 1.0 | Full | 0.6524726851411089 | 0.06777572101415535 | 5 |

Ordering statements:
- 50% completed checkpoint only: Clinical-only − Full = `0.04499323410013534`; ordering changed versus stored 100% = `no`, but this is not the pre-specified completed estimand.
- 75%: ordering change = `unavailable`; no fit was started.
- Stored 100% Table 7: Clinical-only − Full = `0.10734816251880197`; 95% CI = `[0.03700289487664033, 0.17268175928129118]`; ordering changed versus itself = `no`.

No learning-curve trend is classified because the reduced design did not complete.
