"""Independent NumPy matrices and CLI checker for run363's deletion tournament."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
sys.path.insert(0,str(REPO))
import hashlib,json,sqlite3,subprocess
import numpy as np
from check_circuit import right_shift
from search13_fulltail_invariant_gauge_rank import circuit_matrix

r=json.loads((ROOT/'deletion_result.json').read_text())
f=json.loads((ROOT/'deletion_frozen.json').read_text())
assert hashlib.sha256((REPO/f['source']).read_bytes()).hexdigest()==f['source_sha256']
planned={x['deleted_slot']:tuple(map(tuple,x['native_schedule'])) for x in f['proposals']}
assert set(planned)==set(range(13))
assert r['short_count']==13 and r['deep_count']==3
rank=[x['deleted_slot'] for x in sorted(r['short_rows'],key=lambda x:(x['loss'],x['deleted_slot']))[:3]]
assert rank==[x['deleted_slot'] for x in r['deep_rows']]
records=[]
for stage,key in [('short','short_rows'),('deep','deep_rows')]:
 for row in r[key]:
  path=ROOT/row['candidate']; c=json.loads(path.read_text())
  actual=tuple((g['control'],g['target']) for g in c['gates'] if g['gate']=='cx')
  assert actual==planned[row['deleted_slot']]==tuple(map(tuple,row['native_schedule']))
  assert len(c['gates'])==68 and sum(g['gate']=='u3' for g in c['gates'])==56
  u=circuit_matrix(c['gates']); d=u@right_shift(4)@u.conj().T
  off=d-np.diag(np.diag(d)); loss=float(np.sum(abs(off)**2)/16)
  p=subprocess.run(['python3',str(REPO/'check_circuit.py'),str(path)],capture_output=True,text=True)
  q=json.loads(p.stdout)
  assert abs(loss-row['loss'])<1e-11
  assert abs(float(abs(off).max())-q['off_diagonal_error'])<1e-11
  assert abs(q['off_diagonal_error']-row['off_diagonal_error'])<1e-11
  assert q['cnot_count']==12 and q['valid_diagonalizer']==row['valid_diagonalizer']
  records.append({'stage':stage,'deleted_slot':row['deleted_slot'],'candidate':row['candidate'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'numpy_loss':loss,'loss_record_error':abs(loss-row['loss']),'checker_exit':p.returncode,'checker':q})
with sqlite3.connect(REPO/'experiments.sqlite3') as db:
 sha=db.execute('select code_sha256 from attempts where id=363').fetchone()[0]
assert hashlib.sha256((ROOT/'deletion_tournament.py').read_bytes()).hexdigest()==sha
out={'run':363,'source_sha256':f['source_sha256'],'ranking':rank,'code_snapshot_matches':True,'records':records,'maximum_loss_record_error':max(x['loss_record_error'] for x in records)}
(ROOT/'deletion_validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'count':len(records),'valid':sum(x['checker']['valid_diagonalizer'] for x in records),'ranking':rank,'maximum_loss_record_error':out['maximum_loss_record_error']}))
