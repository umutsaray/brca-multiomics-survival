#!/usr/bin/env python3
"""Shared leakage-safe helpers for the reviewer-requested Phase-1 analyses."""
from __future__ import annotations
import hashlib, json, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from lifelines import CoxPHFitter
from sklearn.compose import ColumnTransformer
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from impute_and_split import CLINICAL_COVARIATES, load_clinical_covariates, scan_train_methylation_variance
from prepare_dataset import barcode
from select_features import BLOCKS, rank_block
from cross_validate import train_dl_fold

PREFIXES={"rna":"rna_","mirna":"mirna_","methylation":"meth_","clinical":"clin_"}

def read_frame(path): return pq.read_table(path).to_pandas().set_index("sample_key")

def load_cohort(root: Path, cohort: str = 'mixed'):
    yt,ys=read_frame(root/'y_train.parquet'),read_frame(root/'y_test.parquet')
    y=pd.concat([yt,ys]); schema=pq.ParquetFile(root/'X_train.parquet').schema_arrow.names
    cols=['sample_key',*[c for c in schema if c.startswith(('rna_','mirna_'))]]
    x=pd.concat([pq.read_table(root/'X_train.parquet',columns=cols).to_pandas(),pq.read_table(root/'X_test.parquet',columns=cols).to_pandas()]).set_index('sample_key').loc[y.index]
    y=y.loc[x.index]
    if cohort == 'primary':
        payload=json.loads((root/'primary_tumor_sample_keys.json').read_text(encoding='utf-8'))
        keys=pd.Index(payload['sample_keys']);x=x.loc[keys];y=y.loc[keys]
    elif cohort != 'mixed': raise ValueError(f'Unknown cohort: {cohort}')
    return x,y

def grouped_folds(y,n_splits=5,seed=42):
    groups=pd.Index(y.index).str[:12]
    return list(StratifiedGroupKFold(n_splits=n_splits,shuffle=True,random_state=seed).split(np.zeros(len(y)),y.event,groups))

def collapse_stage(v):
    if pd.isna(v): return v
    s=str(v).strip(); u=s.upper().replace('STAGE ','').replace('STAGE','').strip()
    if u in {'I','IA','IB'}: return 'Stage I'
    if u in {'II','IIA','IIB'}: return 'Stage II'
    if u in {'III','IIIA','IIIB','IIIC'}: return 'Stage III'
    if u=='IV': return 'Stage IV'
    return s

def fold_clinical(root,all_keys,train_keys,test_keys):
    raw=load_clinical_covariates(root/'TCGA-BRCA.clinical.tsv.gz',all_keys)
    stage='ajcc_pathologic_stage.diagnoses'
    print('Gerçek AJCC kategorileri:',sorted(raw[stage].dropna().astype(str).unique()))
    raw[stage]=raw[stage].map(collapse_stage)
    numeric=['age_at_index.demographic']; categorical=[c for c in CLINICAL_COVARIATES if c not in numeric]
    proc=ColumnTransformer([('numeric',SimpleImputer(strategy='median'),numeric),('categorical',Pipeline([('imputer',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float32))]),categorical)])
    a=proc.fit_transform(raw.loc[train_keys]); b=proc.transform(raw.loc[test_keys]); names=['clin_'+x for x in proc.get_feature_names_out()]
    return pd.DataFrame(a,index=train_keys,columns=names,dtype='float32'),pd.DataFrame(b,index=test_keys,columns=names,dtype='float32')

def load_selected_methylation(path,selected,all_keys):
    header=pd.read_csv(path,sep='\t',nrows=0).columns.astype(str).tolist(); fid=header[0]
    raw=[c for c in header[1:] if barcode(c) in set(all_keys)]
    probes=[c.removeprefix('meth_') for c in selected]
    df=pd.read_csv(path,sep='\t',usecols=[fid,*raw],low_memory=False)
    df=df[df[fid].astype(str).isin(probes)].set_index(fid).T
    df.index=pd.Index([barcode(x) for x in df.index],name='sample_key')
    if df.index.has_duplicates: df=df.groupby(level=0,sort=False).mean()
    df.columns=['meth_'+str(c) for c in df.columns]
    return df.reindex(all_keys).reindex(columns=selected).astype('float32')

def prepare_fold(root,base,y,train_idx,test_idx,meth_n=10000,knn=5,chunk=5000):
    trk,tek=y.index[train_idx],y.index[test_idx]
    signature=hashlib.sha256(('\n'.join(trk)+'|'+str(meth_n)+'|'+str(knn)).encode()).hexdigest()[:16]
    cache=root/'phase1_fold_cache';cache.mkdir(exist_ok=True);tp=cache/f'{signature}_train.parquet';ep=cache/f'{signature}_test.parquet'
    if tp.is_file() and ep.is_file():
        print(f'Fold preprocessing cache kullanıldı: {signature}')
        return pd.read_parquet(tp),pd.read_parquet(ep),{'cache':True,'signature':signature}
    chosen,total,seconds=scan_train_methylation_variance(root/'TCGA-BRCA.methylation450.tsv.gz',trk,y.index,meth_n,chunk)
    meth=load_selected_methylation(root/'TCGA-BRCA.methylation450.tsv.gz',chosen,y.index)
    imp=KNNImputer(n_neighbors=knn,keep_empty_features=True); mt=imp.fit_transform(meth.loc[trk]); me=imp.transform(meth.loc[tek])
    mt=pd.DataFrame(mt,index=trk,columns=chosen,dtype='float32'); me=pd.DataFrame(me,index=tek,columns=chosen,dtype='float32')
    ct,ce=fold_clinical(root,y.index,trk,tek)
    train=pd.concat([base.loc[trk],mt,ct],axis=1);test=pd.concat([base.loc[tek],me,ce],axis=1);train.to_parquet(tp);test.to_parquet(ep)
    return train,test,{'meth_total':total,'meth_selected':len(chosen),'variance_seconds':seconds,'cache':False,'signature':signature}

def rank_fold(xtr,ytr,chunk=512):
    rankings={}
    for block,prefix in BLOCKS.items():
        cols=[c for c in xtr if c.startswith(prefix)]
        rankings[block]=rank_block(xtr,cols,ytr,chunk,50,1e-7)
    return rankings

def select_from_rankings(rankings,xtr,top_n=200,fdr=True):
    from statsmodels.stats.multitest import multipletests
    selected=[]; stats={}
    for block,r in rankings.items():
        q=r.copy(); valid=q.p_value.notna(); q['fdr_bh']=np.nan
        q.loc[valid,'fdr_bh']=multipletests(q.loc[valid,'p_value'],method='fdr_bh')[1]
        eligible=q[q.fdr_bh<0.05] if fdr else q[valid]
        chosen=eligible.head(min(top_n,len(eligible))).feature.tolist(); selected+=chosen
        stats[block]={'tested':int(valid.sum()),'fdr_lt_0_05':int((q.fdr_bh<0.05).sum()),'selected':len(chosen)}
    clinical=[c for c in xtr if c.startswith('clin_')]; selected+=clinical
    return selected,stats

def fit_models(xtr,ytr,xte,yte,args,seed):
    out={}; event=ytr.event.astype(bool); ev=xtr.loc[event].var(); ce=xtr.loc[~event].var(); cols=[c for c in xtr if ev[c]>1e-12 and ce[c]>1e-12]
    sc=StandardScaler(); a=sc.fit_transform(xtr[cols]); b=sc.transform(xte[cols]); d=pd.DataFrame(a,columns=cols); d['time']=ytr.time.to_numpy(); d['event']=ytr.event.to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always'); cox=CoxPHFitter(penalizer=args.cox_penalizer); cox.fit(d,duration_col='time',event_col='event',fit_options={'max_steps':args.cox_max_steps})
    cr=cox.predict_partial_hazard(pd.DataFrame(b,columns=cols)).to_numpy().ravel(); out['coxph']={'risk':cr,'cindex':float(concordance_index_censored(yte.event.astype(bool),yte.time,cr)[0]),'warnings':[str(w.message) for w in caught]}
    rsf=RandomSurvivalForest(n_estimators=args.rsf_trees,min_samples_split=10,min_samples_leaf=5,random_state=seed,n_jobs=-1); rsf.fit(xtr.to_numpy('float32'),Surv.from_arrays(ytr.event.astype(bool),ytr.time)); rr=rsf.predict(xte.to_numpy('float32')); out['rsf']={'risk':rr,'cindex':float(concordance_index_censored(yte.event.astype(bool),yte.time,rr)[0]),'model':rsf}
    score,epoch,val,risk=train_dl_fold_with_risk(xtr,ytr,xte,yte,args,seed); out['deepsurv']={'risk':risk,'cindex':score,'best_epoch':epoch,'validation_cindex':val}
    return out

def train_dl_fold_with_risk(xtr,ytr,xte,yte,args,seed):
    # Reuse the validated trainer, then deterministically repeat prediction-capable training is not exposed.
    # Phase scripts needing OOF risks use Cox/RSF directly and train_explainable_dl-compatible code below.
    import torch
    from sklearn.model_selection import train_test_split
    from train_explainable_dl import OmicsDeepSurv,make_blocks,predict_risk,cox_negative_log_likelihood,cindex,set_seed
    set_seed(seed); by={b:[c for c in xtr if c.startswith(p)] for b,p in PREFIXES.items()}; order=[c for b in PREFIXES for c in by[b]]
    slices={}; cur=0
    for b in PREFIXES: slices[b]=(cur,cur+len(by[b])); cur+=len(by[b])
    fi,vi=train_test_split(np.arange(len(xtr)),test_size=args.validation_size,random_state=seed,stratify=ytr.event)
    A=xtr[order].to_numpy('float32'); B=xte[order].to_numpy('float32'); mu=A[fi].mean(0); sd=A[fi].std(0); sd[sd==0]=1; A=(A-mu)/sd; B=(B-mu)/sd
    fb=make_blocks(A[fi],slices,torch.device('cpu')); vb=make_blocks(A[vi],slices,torch.device('cpu')); tb=make_blocks(B,slices,torch.device('cpu'))
    model=OmicsDeepSurv({b:len(by[b]) for b in PREFIXES},dropout=args.dl_dropout); opt=torch.optim.AdamW(model.parameters(),lr=args.dl_learning_rate,weight_decay=args.dl_weight_decay)
    t=torch.tensor(ytr.iloc[fi].time.to_numpy(),dtype=torch.float32); e=torch.tensor(ytr.iloc[fi].event.to_numpy(),dtype=torch.float32); best=None; bs=-np.inf; be=0; wait=0
    for epoch in range(1,args.dl_epochs+1):
        model.train(); opt.zero_grad(); loss=cox_negative_log_likelihood(model(fb),t,e); loss.backward(); opt.step(); vr=predict_risk(model,vb); s=cindex(ytr.iloc[vi].event.to_numpy(),ytr.iloc[vi].time.to_numpy(),vr)
        if s>bs+1e-6: bs=s;be=epoch;wait=0;best={k:v.detach().clone() for k,v in model.state_dict().items()}
        else: wait+=1
        if wait>=args.dl_patience: break
    model.load_state_dict(best); risk=predict_risk(model,tb); return float(concordance_index_censored(yte.event.astype(bool),yte.time,risk)[0]),be,bs,risk

def add_model_args(p):
    p.add_argument('--cox-penalizer',type=float,default=.1);p.add_argument('--cox-max-steps',type=int,default=50);p.add_argument('--rsf-trees',type=int,default=500);p.add_argument('--dl-epochs',type=int,default=300);p.add_argument('--dl-patience',type=int,default=30);p.add_argument('--dl-learning-rate',type=float,default=1e-3);p.add_argument('--dl-weight-decay',type=float,default=1e-4);p.add_argument('--dl-dropout',type=float,default=.3);p.add_argument('--validation-size',type=float,default=.15)
