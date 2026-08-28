# Additional Metrics v2

sksurv cumulative_dynamic_auc survival-function-weighted mean over the supplied 1/3/5-year grid; not arithmetic mean.
IBS grid: 100 equally spaced points from day 365 through day 1825.
IPCW censoring source: for each held-out outer fold, the other four OOF-fold outcomes; never that held fold.

## mixed

| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |
|---|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | coxph | 0.6279095216996384 | nan | 0.724780869355947 | 0.7106857810987027 | nan | 0.08881825609898397 |
| Clinical-only | rsf | 0.6294629356333962 | nan | 0.7328130115830087 | 0.6906787894697782 | nan | 0.08030998730826099 |
| Clinical-only | deepsurv | 0.6310334831205384 | nan | 0.6445688145946818 | 0.5840129608155799 | nan | 0.09163552081692851 |
| Full | coxph | 0.6471561898727105 | nan | 0.7172701040309727 | 0.6283336853886183 | nan | 0.09224838460573087 |
| Full | rsf | 0.6563157987996423 | nan | 0.6887239374252337 | 0.6592146764275363 | nan | 0.08893523588573252 |
| Full | deepsurv | 0.6089233216566455 | nan | 0.6485252858299404 | 0.6652065711758844 | nan | 0.09385635052785903 |

Five-year calibration group details (including Greenwood CI, unique participants at risk, and events) are preserved without rounding in `additional_metrics_v2.json`.

## primary

| Configuration | Model | Uno C | AUC 1y | AUC 3y | AUC 5y | Weighted mean AUC | IBS |
|---|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | coxph | 0.5309700680307541 | 0.891449904539147 | 0.7769179575205037 | 0.7322205313472967 | 0.7613720734344559 | 0.08125873432055347 |
| Clinical-only | rsf | 0.5324908848611061 | 0.8838506384144744 | 0.766565404936645 | 0.7167684902510854 | 0.7496644108065019 | 0.07363326200187073 |
| Clinical-only | deepsurv | 0.5271477709467122 | 0.9119971256500268 | 0.7519963852756352 | 0.7126325648039296 | 0.7448520864144466 | 0.07683708066956368 |
| Full | coxph | 0.5147007286569559 | 0.7342009491105498 | 0.6623165254547048 | 0.6278452283008362 | 0.6511415412546482 | 0.08326204245642184 |
| Full | rsf | 0.5134777183022649 | 0.6489748499807766 | 0.6649823851795074 | 0.6291417882176457 | 0.6254106920133089 | 0.08267871845159261 |
| Full | deepsurv | 0.4270524349475583 | 0.5719288775198998 | 0.5550476462054192 | 0.5137756039112213 | 0.5284399672863609 | 0.08341696250240421 |

Five-year calibration group details (including Greenwood CI, unique participants at risk, and events) are preserved without rounding in `additional_metrics_v2.json`.
