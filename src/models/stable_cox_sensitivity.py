#!/usr/bin/env python3
"""Clinical-only elastic-net Cox sensitivity after rank-deficiency diagnosis."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,pyarrow.parquet as pq
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from phase1_common import load_cohort,grouped_folds
def main():
 root=Path(__file__).resolve().parent;out={}
 for cohort in ('mixed','primary'):
  _,y=load_cohort(root,cohort);scores=[]
  for fold,(tr,te) in enumerate(grouped_folds(y,5,42),1):
   sig=hashlib.sha256(('\n'.join(y.index[tr])+'|10000|5').encode()).hexdigest()[:16];tp=root/'phase1_fold_cache'/f'{sig}_train.parquet';ep=root/'phase1_fold_cache'/f'{sig}_test.parquet';names=pq.ParquetFile(tp).schema_arrow.names;cols=[c for c in names if c.startswith('clin_')];a=pq.read_table(tp,columns=cols).to_pandas();b=pq.read_table(ep,columns=cols).to_pandas();sc=StandardScaler();u=sc.fit_transform(a);v=sc.transform(b);model=CoxnetSurvivalAnalysis(l1_ratio=.5,alphas=[.01],max_iter=100000,tol=1e-7);model.fit(u,Surv.from_arrays(y.iloc[tr].event.astype(bool),y.iloc[tr].time));risk=model.predict(v);score=float(concordance_index_censored(y.iloc[te].event.astype(bool),y.iloc[te].time,risk)[0]);scores.append(score);print(cohort,fold,score)
  out[cohort]={'l1_ratio':.5,'alpha':.01,'fold_cindices':scores,'mean':float(np.mean(scores)),'sample_sd':float(np.std(scores,ddof=1))}
 (root/'stable_cox_sensitivity.json').write_text(json.dumps(out,indent=2),encoding='utf-8');(root/'stable_cox_sensitivity.md').write_text('# Stable Cox Sensitivity\n\n```json\n'+json.dumps(out,indent=2)+'\n```\n',encoding='utf-8')
if __name__=='__main__':main()
