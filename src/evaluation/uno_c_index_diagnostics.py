#!/usr/bin/env python3
"""Diagnose fold-local censoring weights used by Uno's C-index."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sksurv.nonparametric import CensoringDistributionEstimator
from sksurv.util import Surv

MODELS = ("coxph", "rsf", "deepsurv")

def markdown_table(frame: pd.DataFrame) -> str:
    cols=list(frame.columns);lines=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
    for row in frame.itertuples(index=False,name=None):lines.append('| '+' | '.join(str(x) for x in row)+' |')
    return '\n'.join(lines)

def main() -> None:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--output',type=Path,default=Path('uno_c_index_diagnostics.md'));a=p.parse_args()
    rows=[]
    for cohort,filename in (("mixed","ablation_oof_extended.parquet"),("primary","ablation_oof_primary.parquet")):
        data=pd.read_parquet(a.input_dir/filename)
        for config in ("Clinical-only","Full"):
            d=data[data.configuration.eq(config)]
            for fold in sorted(d.fold.unique()):
                train=d[d.fold.ne(fold)];test=d[d.fold.eq(fold)];sy=Surv.from_arrays(train.event.astype(bool),train.time)
                estimator=CensoringDistributionEstimator().fit(sy)
                ipcw=estimator.predict_ipcw(Surv.from_arrays(test.event.astype(bool),test.time))
                positive=ipcw[ipcw>0];tau=float(test.loc[test.event.eq(1),'time'].max())
                rows.append({'cohort':cohort,'configuration':config,'fold':int(fold),'n_test':len(test),'effective_max_event_time_days':tau,'max_observed_time_days':float(test.time.max()),'at_risk_at_effective_horizon':int((test.time>=tau).sum()),'events_by_effective_horizon':int(((test.event==1)&(test.time<=tau)).sum()),'event_ipcw_min':float(positive.min()) if len(positive) else np.nan,'event_ipcw_median':float(np.median(positive)) if len(positive) else np.nan,'event_ipcw_max':float(positive.max()) if len(positive) else np.nan,'concordance_weight_max_ipcw_squared':float(np.square(positive).max()) if len(positive) else np.nan,'n_nonzero_event_weights':len(positive)})
    out=pd.DataFrame(rows)
    lines=['# Uno C-index diagnostics','', '- Script: `src/evaluation/uno_c_index_diagnostics.py`', '- Inputs: retained outer-fold OOF predictions from `ablation_oof_extended.parquet` and `ablation_oof_primary.parquet`.', '- Fold definition: the stored five participant-grouped outer folds (global split seed 42).', '- The retained analysis called `concordance_index_ipcw` with `tau=None`, so no fixed truncation horizon was imposed. `effective_max_event_time_days` is the latest held-fold event time that can contribute; maximum observed follow-up is reported separately.', '- IPCW values are estimated from the other four folds; zero entries correspond to censored observations. The concordance implementation applies squared IPCW event weights.', '',markdown_table(out),'']
    grouped=out.groupby(['cohort','configuration']).agg(max_event_ipcw=('event_ipcw_max','max'),max_squared_weight=('concordance_weight_max_ipcw_squared','max')).reset_index()
    lines += ['## Observed weight ranges','',markdown_table(grouped),'','Large weights are reported descriptively only. Association with the Harrell–Uno discrepancy is not asserted without a dedicated perturbation analysis.','']
    a.output.write_text('\n'.join(lines),encoding='utf-8');print(out.to_string(index=False));print(f'Wrote {a.output}')
if __name__=='__main__':main()
