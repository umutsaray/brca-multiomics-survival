#!/usr/bin/env python3
"""Additional survival metrics from leakage-safe out-of-fold predictions."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
from sksurv.metrics import concordance_index_ipcw,cumulative_dynamic_auc,integrated_brier_score
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--times-days',type=float,nargs='+',default=[365.,1095.,1825.]);a=p.parse_args();root=a.input_dir.resolve();d=pd.read_parquet(root/'phase1_oof_predictions.parquet');surv=Surv.from_arrays(d.event.astype(bool),d.time);times=np.asarray(a.times_days);result={}
 for m in ['coxph','rsf','deepsurv']:
  uno=float(concordance_index_ipcw(surv,surv,d[m].to_numpy())[0]);aucs,mean_auc=cumulative_dynamic_auc(surv,surv,d[m].to_numpy(),times);result[m]={'uno_c_index':uno,'time_dependent_auc':dict(zip(map(str,times),map(float,aucs))),'mean_auc':float(mean_auc)};print(m,result[m])
 # Cross-fit a one-covariate Cox calibration layer: each held-out fold is calibrated
 # only with OOF risks/outcomes from the other four folds.
 grid=np.linspace(times.min(),times.max(),100); estimates={m:np.full((len(d),len(grid)),np.nan) for m in result}
 for held in sorted(d.fold.unique()):
  fit=d.fold.ne(held);test=d.fold.eq(held)
  for m in result:
   z=d.loc[fit,[m,'time','event']].rename(columns={m:'risk'});cal=CoxPHFitter(penalizer=.1);cal.fit(z,'time','event');sf=cal.predict_survival_function(d.loc[test,[m]].rename(columns={m:'risk'}),times=grid).T.to_numpy();estimates[m][np.flatnonzero(test),:]=sf
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(7,5))
 for m,est in estimates.items():
  result[m]['integrated_brier_score']=float(integrated_brier_score(surv,surv,est,grid));pred=est[:,-1];work=d.assign(pred5=pred,bin=pd.qcut(pred,5,duplicates='drop'));pts=[]
  for _,g in work.groupby('bin',observed=True):
   kt,ks=kaplan_meier_estimator(g.event.astype(bool),g.time);pos=np.searchsorted(kt,1825.,side='right')-1;obs=float(ks[pos]) if pos>=0 else 1.;pts.append((float(g.pred5.mean()),obs))
  result[m]['calibration_5year_points']=[{'predicted':x,'observed_km':z} for x,z in pts];ax.plot([x for x,_ in pts],[z for _,z in pts],marker='o',label=m)
 ax.plot([0,1],[0,1],'k--');ax.set(xlabel='Predicted 5-year survival',ylabel='Kaplan–Meier observed 5-year survival');ax.legend();fig.tight_layout();fig.savefig(root/'calibration_5year.png',dpi=600);plt.close(fig)
 result['absolute_survival_method']='One-covariate Cox calibration cross-fitted by outer fold; held-fold outcomes were not used to calibrate their predictions.'
 (root/'additional_metrics.json').write_text(json.dumps(result,indent=2),encoding='utf-8');(root/'additional_metrics.md').write_text('# Additional OOF Metrics\n\n```json\n'+json.dumps(result,indent=2)+'\n```\n',encoding='utf-8')
if __name__=='__main__':main()
