#!/usr/bin/env python3
"""Plot retained five-year calibration groups without recomputing statistics."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

COHORTS=(("mixed","Mixed cohort"),("primary","Primary-tumor-only cohort"))
MODELS=(("coxph","CoxPH"),("rsf","RSF"),("deepsurv","DeepSurv"))
CONFIGS=(("Clinical-only","Clinical-only","#0072B2","o"),("Full","Full multi-omics","#D55E00","s"))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dpi',type=int,default=300);a=p.parse_args();data=json.loads(a.input.read_text(encoding='utf-8'))
 fig,axes=plt.subplots(2,3,figsize=(13.2,8.7),sharex=True,sharey=True)
 for row,(cohort,cohort_label) in enumerate(COHORTS):
  for col,(model,model_label) in enumerate(MODELS):
   ax=axes[row,col];ax.plot([0,1],[0,1],ls='--',lw=1.1,color='#666666',label='Perfect calibration')
   for config,label,color,marker in CONFIGS:
    groups=data[cohort][config][model]['calibration_5year_groups']
    valid=[g for g in groups if np.isfinite(g['predicted_survival_mean']) and np.isfinite(g['observed_km_survival']) and all(np.isfinite(g['greenwood_95_ci']))]
    if not valid:
     ax.text(.5,.5,f'{label}: not available',ha='center',va='center',transform=ax.transAxes,color=color);continue
    x=np.array([g['predicted_survival_mean'] for g in valid]);y=np.array([g['observed_km_survival'] for g in valid]);lo=np.array([g['greenwood_95_ci'][0] for g in valid]);hi=np.array([g['greenwood_95_ci'][1] for g in valid]);order=np.argsort(x);x,y,lo,hi=x[order],y[order],lo[order],hi[order];ordered=[valid[i] for i in order]
    ax.errorbar(x,y,yerr=np.vstack((y-lo,hi-y)),fmt=marker+'-',color=color,lw=1.4,ms=5,capsize=2.5,label=label)
    for xx,yy,g in zip(x,y,ordered):ax.annotate(f"r={g['at_risk_unique_participants_at_5y']}",(xx,yy),xytext=(4,4),textcoords='offset points',fontsize=6.5,color=color)
   ax.set_title(f'{cohort_label} — {model_label}',fontsize=10,fontweight='bold');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.grid(alpha=.2,lw=.6)
   if row==1:ax.set_xlabel('Predicted 5-year survival probability')
   if col==0:ax.set_ylabel('Observed 5-year KM survival')
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,.985))
 fig.suptitle('Five-year calibration from participant-grouped out-of-fold predictions',fontsize=13,y=.998)
 fig.text(.5,.012,'Error bars: Greenwood 95% CI. Point labels r = unique participants at risk at 5 years. Mixed Full CoxPH has four retained groups because risk-group binning produced four distinct bins.',ha='center',fontsize=8)
 fig.tight_layout(rect=(0,.04,1,.94));a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=a.dpi,bbox_inches='tight',facecolor='white');plt.close(fig);print(f'Wrote {a.output} ({a.dpi} DPI; 6 panels, 2 configurations per panel)')
if __name__=='__main__':main()
