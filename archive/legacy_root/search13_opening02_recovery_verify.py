"""Independent saved-gate and CLI validation of the run354 recovery."""
import hashlib
import json
import sqlite3
import subprocess
from pathlib import Path
import numpy as np
from check_circuit import right_shift
from search13_fulltail_invariant_gauge_rank import circuit_matrix

root=Path(__file__).resolve().parent
r=json.loads((root/'search13_opening02_recovery_result.json').read_text())
f=json.loads((root/r['frozen_proposals']).read_text())
planned={x['proposal']:tuple(map(tuple,x['schedule'])) for x in f['proposals']}
assert len(planned)==60 and len(f['omitted'])==5
assert r['short_count']==60 and r['deep_count']==6
ranking=[x['proposal'] for x in sorted(r['short_rows'],key=lambda x:(x['loss'],x['proposal']))[:6]]
assert ranking==[x['proposal'] for x in r['deep_rows']]
assert sum('recovered_from_run' in x for x in r['short_rows'])==7
records=[]
for stage,key in [('short','short_rows'),('deep','deep_rows')]:
 for x in r[key]:
  path=root/x['candidate']; cand=json.loads(path.read_text())
  schedule=tuple((g['control'],g['target']) for g in cand['gates'] if g['gate']=='cx')
  assert schedule==planned[x['proposal']]==tuple(map(tuple,x['schedule']))
  p=subprocess.run(['python3',str(root/'check_circuit.py'),str(path)],capture_output=True,text=True)
  q=json.loads(p.stdout)
  u=circuit_matrix(cand['gates']); d=u@right_shift(4)@u.conj().T
  loss=float(np.sum(np.abs(d-np.diag(np.diag(d)))**2)/16)
  assert abs(loss-x['loss'])<1e-11
  assert abs(q['off_diagonal_error']-x['off_diagonal_error'])<1e-11
  assert q['cnot_count']==13
  if 'candidate_sha256' in x: assert hashlib.sha256(path.read_bytes()).hexdigest()==x['candidate_sha256']
  records.append({'stage':stage,'proposal':x['proposal'],'candidate':x['candidate'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'numpy_loss':loss,'loss_record_error':abs(loss-x['loss']),'cli_exit':p.returncode,'checker':q})
with sqlite3.connect(root/'experiments.sqlite3') as db:
 sha=db.execute('select code_sha256 from attempts where id=354').fetchone()[0]
assert hashlib.sha256((root/'search13_opening02_recovery.py').read_bytes()).hexdigest()==sha
out={'run':354,'source_sha256':r['source_sha256'],'ranking':ranking,'recovered_count':7,'new_short_count':53,'deep_count':6,'code_snapshot_matches':True,'records':records,'maximum_loss_record_error':max(x['loss_record_error'] for x in records)}
(root/'search13_opening02_recovery_validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'count':len(records),'valid':sum(x['checker']['valid_diagonalizer'] for x in records),'ranking':ranking,'maximum_loss_record_error':out['maximum_loss_record_error']}))
