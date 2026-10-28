#!/usr/bin/env python3
"""Detached orchestration for resumable Phase 17 Task 6, reporting, commit, and push."""
from __future__ import annotations
import json,subprocess,sys,time
from datetime import datetime
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'results'/'phase17'
ANALYSIS=[sys.executable,str(REPO/'src/evaluation/phase17_task6_learning_curve.py'),'--input-dir',str(REPO.parent),'--repo-dir',str(REPO),'--output-dir',str(OUT),'--outer-seed','42','--meth-prefilter','10000','--methylation-chunksize','5000','--knn-neighbors','5','--top-n','200','--cox-penalizer','0.1','--cox-max-steps','50','--bootstrap-iterations','9999','--bootstrap-seed','42','--max-projected-remaining-seconds','25200','--finalization-reserve-seconds','600']
SUMMARY=[sys.executable,str(REPO/'src/evaluation/phase17_summary.py'),'--results-dir',str(OUT)]

def run(command):
 print('RUN',subprocess.list2cmdline(command),flush=True);return subprocess.run(command,cwd=REPO,check=False).returncode

def main():
 started=datetime.now().astimezone().isoformat();rc=run(ANALYSIS)
 if rc!=0:print(f'ANALYSIS_FAILED returncode={rc}',flush=True);raise SystemExit(rc)
 payload=json.loads((OUT/'task6_learning_curve.json').read_text(encoding='utf-8'));status=payload['status']
 rc=run(SUMMARY)
 if rc!=0:print(f'SUMMARY_FAILED returncode={rc}',flush=True);raise SystemExit(rc)
 background={'started_at':started,'finished_at':datetime.now().astimezone().isoformat(),'analysis_returncode':0,'analysis_status':status,'analysis_command':subprocess.list2cmdline(ANALYSIS),'summary_command':subprocess.list2cmdline(SUMMARY),'commit_and_push_requested':True}
 (OUT/'task6_background_run.json').write_text(json.dumps(background,indent=2),encoding='utf-8')
 paths=['src/evaluation/phase17_task6_learning_curve.py','src/evaluation/phase17_task6_background_runner.py','results/phase17/task6_checkpoints','results/phase17/task6_learning_curve.json','results/phase17/task6_learning_curve.md','results/phase17/task6_progress.json','results/phase17/task6_background_run.json','results/phase17/phase17_summary.md']
 rc=run(['git','add','--',*paths])
 if rc!=0:raise SystemExit(rc)
 message='Complete resumed Phase 17 Task 6 learning-curve analysis' if status=='completed' else 'Checkpoint Phase 17 Task 6 after 7-hour runtime stop'
 changed=subprocess.run(['git','diff','--cached','--quiet'],cwd=REPO,check=False).returncode
 if changed:
  rc=run(['git','commit','-m',message])
  if rc!=0:raise SystemExit(rc)
 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
 for attempt in range(1,11):
  rc=run(['git','push','origin','HEAD'])
  if rc==0:break
  print(f'PUSH_RETRY attempt={attempt} commit={commit}',flush=True);time.sleep(30)
 if rc!=0:print(f'PUSH_FAILED commit={commit}',flush=True);raise SystemExit(rc)
 print(f'PHASE17_TASK6_STATUS={status}',flush=True);print(f'PHASE17_TASK6_COMMIT={commit}',flush=True)

if __name__=='__main__':main()
