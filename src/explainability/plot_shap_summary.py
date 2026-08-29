#!/usr/bin/env python3
"""Plot the retained corrected-cohort aggregate SHAP top-20 ranking."""
from __future__ import annotations
import argparse,json,textwrap
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dpi',type=int,default=300);a=p.parse_args();data=json.loads(a.input.read_text(encoding='utf-8'));rows=data['aggregate_top20']
 if len(rows)!=20:raise ValueError(f'Expected 20 aggregate features, found {len(rows)}')
 features=[r['feature'] for r in rows][::-1];values=np.array([r['mean_abs_shap_across_folds_where_selected'] for r in rows][::-1]);folds=np.array([r['folds_selected'] for r in rows][::-1]);norm=mpl.colors.Normalize(vmin=1,vmax=5);cmap=mpl.colormaps['viridis'];colors=cmap(norm(folds))
 fig,ax=plt.subplots(figsize=(11.5,8.5));bars=ax.barh(np.arange(len(rows)),values,color=colors,height=.72)
 ax.set_yticks(np.arange(len(rows)),features,fontsize=8);ax.set_xlabel('Mean |SHAP| across folds where selected');ax.set_ylabel('Feature');ax.grid(axis='x',alpha=.22,lw=.6);ax.set_axisbelow(True)
 offset=max(values)*.012
 for bar,n in zip(bars,folds):ax.text(bar.get_width()+offset,bar.get_y()+bar.get_height()/2,f'{n}/5 folds',va='center',fontsize=7.5)
 ax.set_xlim(0,max(values)*1.19);sm=mpl.cm.ScalarMappable(norm=norm,cmap=cmap);cbar=fig.colorbar(sm,ax=ax,pad=.015,fraction=.025,ticks=[1,2,3,4,5]);cbar.set_label('Folds in whose top 20',fontsize=8)
 ax.set_title('Corrected primary-tumor-only full multi-omics DeepSurv: aggregate SHAP ranking',fontsize=12,fontweight='bold',pad=12)
 jaccard=data['mean_pairwise_jaccard']
 caption=('Corrected primary-tumor-only cohort, full multi-omics DeepSurv model. Feature importance rankings are unstable across cross-validation folds '
          f'(mean pairwise top-20 Jaccard index = {jaccard:.3f}); this figure reflects one representative aggregate ranking, not a stable or reproducible signal.')
 fig.text(.08,.015,'\n'.join(textwrap.wrap(caption,145)),ha='left',va='bottom',fontsize=8)
 fig.tight_layout(rect=(0,.075,1,1));a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=a.dpi,bbox_inches='tight',facecolor='white');plt.close(fig);print(f'Wrote {a.output} ({a.dpi} DPI; {len(rows)} bars; stored mean Jaccard={data["mean_pairwise_jaccard"]})')
if __name__=='__main__':main()
