# OOF Risk-scale Diagnostics

The completed OOF files do not store training-fold risk means. Therefore a leakage-safe training-mean-centered pooled C-index cannot be reconstructed and is not fabricated here. The supplied centered value subtracts each held-out fold's own mean and is labelled diagnostic-only; it must not be used as a performance estimate. Future model fits must persist the training-risk mean and subtract that constant from their held-out risks.

## mixed

| Model | Fold | N | Mean | Sample variance |
|---|---:|---:|---:|---:|
| coxph | 1 | 168 | 33.294518766606494 | 71847.08816545026 |
| coxph | 2 | 168 | 4613.674777521175 | 3569982159.4662733 |
| coxph | 3 | 168 | 33.38856128333034 | 83177.63495792162 |
| coxph | 4 | 169 | 2.2229510140477426 | 17.74381053762658 |
| coxph | 5 | 169 | 1060.0631022713876 | 188879080.67330298 |
| rsf | 1 | 168 | 18.878121104593692 | 144.71961681794247 |
| rsf | 2 | 168 | 19.1031348141559 | 90.18412510778963 |
| rsf | 3 | 168 | 16.61962566473186 | 90.6477610846832 |
| rsf | 4 | 169 | 16.770001306080196 | 89.4499125571923 |
| rsf | 5 | 169 | 17.408051081612214 | 96.02688040592723 |
| deepsurv | 1 | 168 | 0.0795726281856852 | 0.010112180967788143 |
| deepsurv | 2 | 168 | 0.6436978321345079 | 1.4356497552266014 |
| deepsurv | 3 | 168 | -0.15986660757057725 | 1.8388692675597191 |
| deepsurv | 4 | 169 | 0.8431728132642232 | 0.8246122056016001 |
| deepsurv | 5 | 169 | -0.2633692050606189 | 0.2821839222578079 |

| Model | Risk definition | Pooled uncentered C | Test-fold-mean centered C (diagnostic only) |
|---|---|---:|---:|
| coxph | lifelines predict_partial_hazard (partial hazard) | 0.6411477013368112 | 0.5400717313335507 |
| rsf | scikit-survival RandomSurvivalForest.predict risk score | 0.6788392566025432 | 0.6806868818606673 |
| deepsurv | raw neural-network scalar output (log-risk) | 0.6432344310401044 | 0.6411259645690686 |

## primary

| Model | Fold | N | Mean | Sample variance |
|---|---:|---:|---:|---:|
| coxph | 1 | 153 | 641.706865326109 | 60345017.061176516 |
| coxph | 2 | 153 | 225.3313262508023 | 6774531.053631741 |
| coxph | 3 | 152 | 44.326056341609714 | 148949.90214295665 |
| coxph | 4 | 152 | 1288281.9647163074 | 136686335022564.53 |
| coxph | 5 | 152 | 5.458938093022604 | 370.4300443323084 |
| rsf | 1 | 153 | 15.669331805719771 | 73.0803879131428 |
| rsf | 2 | 153 | 14.154510447877138 | 78.27849664155889 |
| rsf | 3 | 152 | 15.20058286724704 | 67.41658640510848 |
| rsf | 4 | 152 | 17.681277376661402 | 78.379668342092 |
| rsf | 5 | 152 | 14.596619836508886 | 78.42817753741794 |
| deepsurv | 1 | 153 | 0.19529121417820064 | 0.007167748667887672 |
| deepsurv | 2 | 153 | -0.11443022653168324 | 0.05866229293074695 |
| deepsurv | 3 | 152 | 0.04220593710871119 | 0.8105765344331961 |
| deepsurv | 4 | 152 | 1.4709474506640905 | 5.397126694453295 |
| deepsurv | 5 | 152 | -0.21509678573592714 | 0.3550873731000774 |

| Model | Risk definition | Pooled uncentered C | Test-fold-mean centered C (diagnostic only) |
|---|---|---:|---:|
| coxph | lifelines predict_partial_hazard (partial hazard) | 0.6713727661406909 | 0.5509519854818369 |
| rsf | scikit-survival RandomSurvivalForest.predict risk score | 0.6495955215157947 | 0.6467349512472702 |
| deepsurv | raw neural-network scalar output (log-risk) | 0.5766663590784658 | 0.5747593122327828 |
