"""Read-only architecture source/candidate/axis audit."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));from verify import matrix,cycle
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'curvature'));from directional_check import axis_matrix
sys.path.insert(0,str(ROOT));from delete14_search import u3_to_axis
folder=ROOT/'collaboration/work/numerical_label12/architecture_corrected';cfg=json.loads((folder/'config.json').read_text());result=json.loads((folder/'result.json').read_text())
source=ROOT/cfg['source'];data=json.loads(source.read_text());u,_=matrix(data);cursor=0;locals_=[]
for slot in range(14):
 for q in range(4):
  g=data['gates'][cursor];cursor+=1;assert g['gate']=='u3' and g['qubit']==q;locals_.append(g)
 if slot<13 and cfg['optimizer_schedule'][slot] is not None:
  c,t=cfg['optimizer_schedule'][slot];assert data['gates'][cursor]=={'gate':'cx','control':c,'target':t};cursor+=1
assert cursor==len(data['gates']);axis=np.array([u3_to_axis(g) for g in locals_]);v=axis_matrix(axis,cfg['optimizer_schedule']);phase=np.vdot(u,v);phase/=abs(phase)
rows=[]
for index,row in enumerate(result['rows']):
 proposal=cfg['proposals'][index];slot=proposal['changed_optimizer_slot'];wanted=cfg['optimizer_schedule'].copy();wanted[slot]=proposal['edge'];assert wanted==proposal['schedule']
 path=ROOT/row['candidate'];d=json.loads(path.read_text());cx=[(g['control'],g['target']) for g in d['gates'] if g['gate']=='cx'];assert cx==[tuple(e) for e in wanted if e is not None]
 m,count=matrix(d);a=axis_matrix(row['axis_checkpoint'],wanted);p=np.vdot(m,a);p/=abs(p);error=float(np.max(abs(a-p*m)));assert error<1e-12
 b=m@cycle()@m.conj().T;off=b-np.diag(np.diag(b));rows.append({'candidate':row['candidate'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'count':count,'axis_saved_gate_error':error,'offloss':float(np.sum(abs(off)**2)/16)})
report={'source_hash_matches':hashlib.sha256(source.read_bytes()).hexdigest()==cfg['source_sha256'],'source_actual':cfg['source'],'source_axis_phase_error':float(np.max(abs(v-phase*u))),'corrected_source_binding':'Actual source377base verified independently before fit.','rows':rows}
assert report['source_hash_matches'];Path(__file__).with_name('independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
