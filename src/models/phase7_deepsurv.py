#!/usr/bin/env python3
"""Leakage-corrected nested DeepSurv for Clinical-only and Full OOF predictions."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd,torch
from sklearn.compose import ColumnTransformer
from sklearn.impute import KNNImputer,SimpleImputer
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sksurv.metrics import concordance_index_censored
from phase1_common import load_cohort,grouped_folds,rank_fold,select_from_rankings,load_selected_methylation
from impute_and_split import CLINICAL_COVARIATES,load_clinical_covariates,scan_train_methylation_variance
from phase1_common import collapse_stage
from train_explainable_dl import PREFIXES,OmicsDeepSurv,FlatOmicsWrapper,make_blocks,predict_risk,cox_negative_log_likelihood,set_seed

def inner_split(y,outer_train_idx,seed):
 sub=y.iloc[outer_train_idx];groups=sub.index.str[:12];splitter=StratifiedGroupKFold(n_splits=7,shuffle=True,random_state=seed);fit_rel,val_rel=next(splitter.split(np.zeros(len(sub)),sub.event,groups));return outer_train_idx[fit_rel],outer_train_idx[val_rel]

def clinical(root,all_keys,fit_keys,val_keys,test_keys):
 raw=load_clinical_covariates(root/'TCGA-BRCA.clinical.tsv.gz',all_keys);stage='ajcc_pathologic_stage.diagnoses';raw[stage]=raw[stage].map(collapse_stage);numeric=['age_at_index.demographic'];categorical=[c for c in CLINICAL_COVARIATES if c not in numeric]
 proc=ColumnTransformer([('numeric',SimpleImputer(strategy='median'),numeric),('categorical',Pipeline([('imputer',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float32))]),categorical)]);proc.fit(raw.loc[fit_keys]);names=['clin_'+x for x in proc.get_feature_names_out()]
 def tx(keys):return pd.DataFrame(proc.transform(raw.loc[keys]),index=keys,columns=names,dtype='float32')
 return tx(fit_keys),tx(val_keys),tx(test_keys)

def prepare_inner(root,base,y,fit_idx,val_idx,test_idx,meth_n,knn,chunk):
 fk,vk,tk=y.index[fit_idx],y.index[val_idx],y.index[test_idx];all_keys=y.index;chosen,total,seconds=scan_train_methylation_variance(root/'TCGA-BRCA.methylation450.tsv.gz',fk,all_keys,meth_n,chunk);meth=load_selected_methylation(root/'TCGA-BRCA.methylation450.tsv.gz',chosen,all_keys);imp=KNNImputer(n_neighbors=knn,keep_empty_features=True);mf=imp.fit_transform(meth.loc[fk]);mv=imp.transform(meth.loc[vk]);mt=imp.transform(meth.loc[tk]);mf=pd.DataFrame(mf,index=fk,columns=chosen,dtype='float32');mv=pd.DataFrame(mv,index=vk,columns=chosen,dtype='float32');mt=pd.DataFrame(mt,index=tk,columns=chosen,dtype='float32');cf,cv,ct=clinical(root,all_keys,fk,vk,tk)
 return pd.concat([base.loc[fk],mf,cf],axis=1),pd.concat([base.loc[vk],mv,cv],axis=1),pd.concat([base.loc[tk],mt,ct],axis=1),{'methylation_total':total,'methylation_prefiltered':len(chosen),'variance_seconds':seconds,'knn_neighbors':knn}

def fit_network(xf,yf,xv,yv,xt,yt,seed,epochs,patience,lr,wd,dropout):
 set_seed(seed);by={b:[c for c in xf if c.startswith(p)] for b,p in PREFIXES.items()};order=[c for b in PREFIXES for c in by[b]];slices={};cur=0
 for b in PREFIXES:slices[b]=(cur,cur+len(by[b]));cur+=len(by[b])
 A=xf[order].to_numpy('float32');V=xv[order].to_numpy('float32');T=xt[order].to_numpy('float32');mu=A.mean(0,dtype=np.float64).astype('float32');sd=A.std(0,dtype=np.float64).astype('float32');sd[sd==0]=1;A=(A-mu)/sd;V=(V-mu)/sd;T=(T-mu)/sd;device=torch.device('cpu');fb=make_blocks(A,slices,device);vb=make_blocks(V,slices,device);tb=make_blocks(T,slices,device);model=OmicsDeepSurv({b:len(by[b]) for b in PREFIXES},dropout=dropout);opt=torch.optim.AdamW(model.parameters(),lr=lr,weight_decay=wd);tt=torch.tensor(yf.time.to_numpy(),dtype=torch.float32);ee=torch.tensor(yf.event.to_numpy(),dtype=torch.float32);best=None;best_val=-np.inf;best_epoch=0;wait=0
 for epoch in range(1,epochs+1):
  model.train();opt.zero_grad();loss=cox_negative_log_likelihood(model(fb),tt,ee);loss.backward();opt.step();vr=predict_risk(model,vb);score=float(concordance_index_censored(yv.event.astype(bool),yv.time,vr)[0])
  if score>best_val+1e-6:best_val=score;best_epoch=epoch;wait=0;best={k:v.detach().clone() for k,v in model.state_dict().items()}
  else:wait+=1
  if wait>=patience:break
 model.load_state_dict(best);risk=predict_risk(model,tb);ci=float(concordance_index_censored(yt.event.astype(bool),yt.time,risk)[0]);return {'risk':risk,'cindex':ci,'best_epoch':best_epoch,'validation_cindex':best_val,'model':model,'order':order,'slices':slices,'fit_scaled':A,'test_scaled':T}

def shap_fold(fit_result,seed,bg_n,explain_n):
 import shap
 model=fit_result['model'];wrapper=FlatOmicsWrapper(model,fit_result['slices']).eval();rng=np.random.default_rng(seed);A=fit_result['fit_scaled'];T=fit_result['test_scaled'];bi=rng.choice(len(A),min(bg_n,len(A)),replace=False);ei=rng.choice(len(T),min(explain_n,len(T)),replace=False);vals=shap.GradientExplainer(wrapper,torch.from_numpy(A[bi])).shap_values(torch.from_numpy(T[ei]));vals=vals[0] if isinstance(vals,list) else vals;vals=np.asarray(vals);vals=vals[...,0] if vals.ndim==3 else vals;imp=np.mean(np.abs(vals),axis=0);idx=np.argsort(imp)[::-1];order=fit_result['order'];return {'background_n':len(bi),'explained_n':len(ei),'top20':[order[i] for i in idx[:20]],'importance':{order[i]:float(imp[i]) for i in idx}}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--random-state',type=int,default=42);p.add_argument('--meth-prefilter',type=int,default=10000);p.add_argument('--knn-neighbors',type=int,default=5);p.add_argument('--methylation-chunksize',type=int,default=5000);p.add_argument('--top-n',type=int,default=200);p.add_argument('--epochs',type=int,default=300);p.add_argument('--patience',type=int,default=30);p.add_argument('--learning-rate',type=float,default=.001);p.add_argument('--weight-decay',type=float,default=.0001);p.add_argument('--dropout',type=float,default=.3);p.add_argument('--background',type=int,default=30);p.add_argument('--explain',type=int,default=30);a=p.parse_args();torch.set_num_threads(1);root=a.input_dir.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);payload={'method':{'outer_split':'5-fold StratifiedGroupKFold, participant-grouped/event-stratified','outer_seed':a.random_state,'internal_split':'first held fold of 7-fold StratifiedGroupKFold on outer-training only (~1/7 validation), participant-grouped/event-stratified','internal_seed_schema':'4200 + outer_fold','training_seed_schema':'42 + outer_fold','selection':'internal-train-only univariate Cox, BH-FDR < 0.05, cap 200/block','all learned preprocessing':'fit internal-train only','torch_threads':1},'cohorts':{}};all_oof=[];shap_folds=[];started=time.perf_counter()
 for cohort in ('mixed','primary'):
  base,y=load_cohort(root,cohort);fold_info=[]
  for fold,(otr,ote) in enumerate(grouped_folds(y,5,a.random_state),1):
   iseed=4200+fold;tseed=a.random_state+fold;fit_idx,val_idx=inner_split(y,otr,iseed);yf,yv,yt=y.iloc[fit_idx],y.iloc[val_idx],y.iloc[ote];cache=root/'phase7_inner_cache';cache.mkdir(exist_ok=True);stem=f'{cohort}_fold{fold}_innerseed{iseed}';fp,vp,tp,mp=[cache/f'{stem}_{x}' for x in ('fit.parquet','validation.parquet','test.parquet','meta.json')]
   if fp.is_file() and vp.is_file() and tp.is_file() and mp.is_file():xf=pd.read_parquet(fp);xv=pd.read_parquet(vp);xt=pd.read_parquet(tp);meta=json.loads(mp.read_text(encoding='utf-8'));stats,prep=meta['selection'],meta['preprocessing'];print('Phase 7 inner-fold cache used:',stem)
   else:
    xf,xv,xt,prep=prepare_inner(root,base,y,fit_idx,val_idx,ote,a.meth_prefilter,a.knn_neighbors,a.methylation_chunksize);rank=rank_fold(xf,yf);selected,stats=select_from_rankings(rank,xf,a.top_n,True);keep=[c for c in selected if c in xf];xf,xv,xt=xf[keep],xv[keep],xt[keep];xf.to_parquet(fp);xv.to_parquet(vp);xt.to_parquet(tp);mp.write_text(json.dumps({'selection':stats,'preprocessing':prep},indent=2),encoding='utf-8')
   clinical_cols=[c for c in xf if c.startswith('clin_')];full_cols=list(xf.columns);fi={'fold':fold,'outer_train_n':len(otr),'internal_train_n':len(fit_idx),'internal_validation_n':len(val_idx),'outer_test_n':len(ote),'outer_seed':a.random_state,'internal_split_seed':iseed,'training_seed':tseed,'selection':stats,'preprocessing':prep,'models':{}}
   for config,cols in (('Clinical-only',clinical_cols),('Full',full_cols)):
    result=fit_network(xf[cols],yf,xv[cols],yv,xt[cols],yt,tseed,a.epochs,a.patience,a.learning_rate,a.weight_decay,a.dropout);fi['models'][config]={'cindex':result['cindex'],'best_epoch':result['best_epoch'],'validation_cindex':result['validation_cindex'],'n_features':len(cols)};print(cohort,fold,config,fi['models'][config])
    for j,key in enumerate(yt.index):all_oof.append({'sample_key':key,'cohort':cohort,'fold':fold,'configuration':config,'event':int(yt.iloc[j].event),'time':float(yt.iloc[j].time),'deepsurv':float(result['risk'][j])})
    if cohort=='primary' and config=='Full':sf=shap_fold(result,tseed,a.background,a.explain);sf.update({'fold':fold,'seed':tseed,'selection':stats,'best_epoch':result['best_epoch']});shap_folds.append(sf)
   fold_info.append(fi)
  payload['cohorts'][cohort]={'n_rows':len(y),'n_participants':int(y.index.str[:12].nunique()),'folds':fold_info}
 pd.DataFrame(all_oof).to_parquet(out/'phase7_deepsurv_oof.parquet',index=False);payload['elapsed_seconds']=time.perf_counter()-started;(out/'phase7_deepsurv_results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
 pairs=[]
 for i in range(5):
  for j in range(i+1,5):A=set(shap_folds[i]['top20']);B=set(shap_folds[j]['top20']);pairs.append({'folds':[i+1,j+1],'intersection':len(A&B),'jaccard':len(A&B)/len(A|B)})
 agg={}
 for f in shap_folds:
  for name,val in f['importance'].items():agg.setdefault(name,[]).append(val)
 ranked=sorted(((k,float(np.mean(v)),len(v)) for k,v in agg.items()),key=lambda x:x[1],reverse=True);sp={'cohort':'primary','pipeline':'Phase 7 internal-train-only preprocessing and selection','explainer':'shap.GradientExplainer','global_seed':a.random_state,'folds':shap_folds,'shap_top20_pairwise':pairs,'mean_pairwise_jaccard':float(np.mean([x['jaccard'] for x in pairs])),'aggregate_top20':[{'feature':k,'mean_abs_shap_across_folds_where_selected':v,'folds_selected':n} for k,v,n in ranked[:20]]};(out/'shap_phase7_results.json').write_text(json.dumps(sp,indent=2),encoding='utf-8');print('Wrote Phase 7 OOF/results/SHAP outputs')
if __name__=='__main__':main()
