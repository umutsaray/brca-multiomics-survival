#!/usr/bin/env python3
"""Build one primary-solid-tumor sample key per TCGA participant."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

def read(path): return pq.read_table(path).to_pandas().set_index('sample_key')
def main():
 root=Path(__file__).resolve().parent
 y=pd.concat([read(root/'y_train.parquet'),read(root/'y_test.parquet')])
 table=y.reset_index();table['participant_id']=table.sample_key.str[:12];table['sample_type']=table.sample_key.str[13:15]
 primary=table[table.sample_type.eq('01')].sort_values(['participant_id','sample_key'])
 counts=primary.groupby('participant_id').size();chosen=primary.drop_duplicates('participant_id',keep='first').set_index('sample_key')
 if not chosen.participant_id.is_unique or not chosen.sample_type.eq('01').all(): raise RuntimeError('Primary cohort invariant failed')
 payload={'source_rows':len(table),'source_participants':int(table.participant_id.nunique()),'source_primary_rows':int((table.sample_type=='01').sum()),'primary_participants':len(chosen),'event_1':int(chosen.event.sum()),'censored_0':int((chosen.event==0).sum()),'removed_from_842':len(table)-len(chosen),'difference_from_762_primary_rows':int((table.sample_type=='01').sum()-len(chosen)),'participants_with_multiple_primary_sample_keys':int((counts>1).sum()),'selection_rule':'Keep sample_type 01 only; sort 15-character sample_key lexicographically and keep the first key per 12-character participant. Vials/aliquots were already collapsed at the 15-character key by prepare_dataset.py using their mean.'}
 (root/'primary_tumor_sample_keys.json').write_text(json.dumps({'summary':payload,'sample_keys':chosen.index.tolist()},indent=2),encoding='utf-8')
 chosen[['participant_id','sample_type','event','time']].to_parquet(root/'primary_tumor_labels.parquet')
 lines=['# Primary-tumor-only Cohort','',f"Selection rule: {payload['selection_rule']}",'']+[f'- {k}: {v}' for k,v in payload.items() if k!='selection_rule']
 (root/'primary_cohort_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print(json.dumps(payload,indent=2))
if __name__=='__main__':main()
