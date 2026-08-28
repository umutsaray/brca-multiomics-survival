#!/usr/bin/env python3
"""Five-fold GradientSHAP analysis from the primary-only leakage-safe outer pipeline."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd,torch,shap
from sklearn.model_selection import train_test_split
from phase1_common import load_cohort,grouped_folds,prepare_fold,rank_fold,select_from_rankings,PREFIXES
from train_explainable_dl import OmicsDeepSurv,FlatOmicsWrapper,make_blocks,predict_risk,cox_negative_log_likelihood,set_seed

def block(f):
 return 'rna' if f.startswith('rna_') else 'mirna' if f.startswith('mirna_') else 'methylation' if f.startswith('meth_') else 'clinical'
def fit_explain(xtr,ytr,xte,seed,bg_n,test_n,epochs,patience):
 set_seed(seed);by={b:[c for c in xtr if c.startswith(p)] for b,p in PREFIXES.items()};order=[c for b in PREFIXES for c in by[b]];slices={};cur=0
 for b in PREFIXES:slices[b]=(cur,cur+len(by[b]));cur+=len(by[b])
 fi,vi=train_test_split(np.arange(len(xtr)),test_size=.15,random_state=seed,stratify=ytr.event);A=xtr[order].to_numpy('float32');B=xte[order].to_numpy('float32');mu=A[fi].mean(0);sd=A[fi].std(0);sd[sd==0]=1;A=(A-mu)/sd;B=(B-mu)/sd
 fb=make_blocks(A[fi],slices,torch.device('cpu'));vb=make_blocks(A[vi],slices,torch.device('cpu'));model=OmicsDeepSurv({b:len(by[b]) for b in PREFIXES},dropout=.3);opt=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001);t=torch.tensor(ytr.iloc[fi].time.to_numpy(),dtype=torch.float32);e=torch.tensor(ytr.iloc[fi].event.to_numpy(),dtype=torch.float32);best=None;score=-np.inf;wait=0
 from sksurv.metrics import concordance_index_censored
 for epoch in range(1,epochs+1):
  model.train();opt.zero_grad();loss=cox_negative_log_likelihood(model(fb),t,e);loss.backward();opt.step();vr=predict_risk(model,vb);s=float(concordance_index_censored(ytr.iloc[vi].event.astype(bool),ytr.iloc[vi].time,vr)[0])
  if s>score+1e-6:score=s;best={k:v.detach().clone() for k,v in model.state_dict().items()};wait=0
  else:wait+=1
  if wait>=patience:break
 model.load_state_dict(best);model.eval();wrapper=FlatOmicsWrapper(model,slices).eval();rng=np.random.default_rng(seed);bi=rng.choice(fi,min(bg_n,len(fi)),replace=False);ei=rng.choice(len(B),min(test_n,len(B)),replace=False);expl=shap.GradientExplainer(wrapper,torch.from_numpy(A[bi]));vals=expl.shap_values(torch.from_numpy(B[ei]));vals=vals[0] if isinstance(vals,list) else vals;vals=np.asarray(vals);vals=vals[...,0] if vals.ndim==3 else vals;return order,np.mean(np.abs(vals),axis=0),len(bi),len(ei),epoch
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--background',type=int,default=30);p.add_argument('--explain',type=int,default=30);p.add_argument('--epochs',type=int,default=300);p.add_argument('--patience',type=int,default=30);p.add_argument('--random-state',type=int,default=42);a=p.parse_args();root=a.input_dir.resolve();base,y=load_cohort(root,'primary');folds=[];start=time.perf_counter()
 for fold,(tr,te) in enumerate(grouped_folds(y,5,a.random_state),1):
  xtr,xte,_=prepare_fold(root,base,y,tr,te);ytr=y.iloc[tr];ranks=rank_fold(xtr,ytr);selected,stats=select_from_rankings(ranks,xtr,200,True);cols=[c for c in selected if c in xtr];order,imp,bn,en,last=fit_explain(xtr[cols],ytr,xte[cols],a.random_state+fold,a.background,a.explain,a.epochs,a.patience);idx=np.argsort(imp)[::-1];top=[order[i] for i in idx[:20]];blocks={}
  for b in PREFIXES:
   mask=np.array([block(f)==b for f in order]);blocks[b]={'feature_count':int(mask.sum()),'total_mean_abs_shap':float(imp[mask].sum()),'mean_abs_shap_per_feature':float(imp[mask].mean()) if mask.any() else None}
  folds.append({'fold':fold,'seed':a.random_state+fold,'selection':stats,'background_n':bn,'explained_n':en,'last_epoch':last,'top20':top,'importance':{order[i]:float(imp[i]) for i in idx},'block_importance':blocks});print('fold',fold,'top20',top)
 pairs=[]
 for i in range(5):
  for j in range(i+1,5):
   A=set(folds[i]['top20']);B=set(folds[j]['top20']);pairs.append({'folds':[i+1,j+1],'intersection':len(A&B),'jaccard':len(A&B)/len(A|B)})
 aggregate={}
 for f in folds:
  for name,val in f['importance'].items():aggregate.setdefault(name,[]).append(val)
 agg=sorted(((k,float(np.mean(v)),len(v)) for k,v in aggregate.items()),key=lambda x:x[1],reverse=True);payload={'cohort':'primary','explainer':'shap.GradientExplainer','background_selection':'random without replacement from neural fit subset','explained_selection':'random without replacement from held-out outer fold','global_seed':a.random_state,'folds':folds,'shap_top20_pairwise':pairs,'mean_pairwise_jaccard':float(np.mean([x['jaccard'] for x in pairs])),'aggregate_top20':[{'feature':k,'mean_abs_shap_across_folds_where_selected':v,'folds_selected':n} for k,v,n in agg[:20]],'elapsed_seconds':time.perf_counter()-start}
 (root/'shap_v2_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');lines=['# SHAP v2 — Primary-tumor-only','',f"Explainer: `{payload['explainer']}`. Background and explained cases were random without replacement with logged fold seeds.",f"Mean pairwise top-20 SHAP Jaccard: {payload['mean_pairwise_jaccard']}.",'','Feature-selection Jaccard and SHAP-ranking Jaccard are different quantities; no equivalence is inferred. Block totals and per-feature means are both reported in JSON to expose feature-count effects.','','| Rank | Feature | Mean | Folds selected |','|---:|---|---:|---:|']
 for i,x in enumerate(payload['aggregate_top20'],1):lines.append(f"| {i} | `{x['feature']}` | {x['mean_abs_shap_across_folds_where_selected']} | {x['folds_selected']} |")
 (root/'shap_v2_results.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print('shap_v2_results.md/json written')
if __name__=='__main__':main()
