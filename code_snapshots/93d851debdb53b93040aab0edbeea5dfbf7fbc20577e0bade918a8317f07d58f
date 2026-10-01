"""Fresh all-pair schedules with independent Haar SU(2) local starts.

Bounded search; unordered duplicate guard covers retained search13 native
candidate files only. Each order uses all six pairs, with no adjacent repeat.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import admissible,fit,serialize
ROOT=Path(__file__).resolve().parent
STEM='search13_balanced_haar_tournament'
SEED=35400

def key(schedule):return tuple(tuple(sorted(e)) for e in schedule)
def haar(rng):
 q=rng.normal(size=(14,4,4));q/=np.linalg.norm(q,axis=-1,keepdims=True)
 s=np.linalg.norm(q[...,1:],axis=-1);angle=2*np.arctan2(s,q[...,0])
 return q[...,1:]*(angle/np.maximum(s,1e-15))[...,None]
def save(case,stage,schedule,angles):
 c=serialize([[] for _ in range(14)],schedule,angles);p=ROOT/f'{STEM}_case{case}_{stage}.json';assert not p.exists();p.write_text(json.dumps(c,indent=2)+'\n');q=evaluate(c);assert q['cnot_count']==13
 return {'candidate':p.name,'cnot_count':13,'off_diagonal_error':q['off_diagonal_error'],'valid_diagonalizer':q['valid_diagonalizer']}
def main():
 out=ROOT/f'{STEM}_result.json';plan=ROOT/f'{STEM}_plan.json';assert not out.exists() and not plan.exists()
 seen=set();coverage=[]
 for p in sorted(ROOT.glob('search13*.json')):
  c=json.loads(p.read_text())
  if isinstance(c,dict) and c.get('n')==4 and 'gates' in c:
   edges=[(g['control'],g['target']) for g in c['gates'] if g['gate']=='cx']
   if len(edges)==13:seen.add(key(edges));coverage.append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 rng=np.random.default_rng(SEED);pairs=[(a,b) for a in range(4) for b in range(a+1,4)];schedules=[]
 for attempt in range(10000):
  if len(schedules)==12:break
  order=[pairs[int(i)] for i in rng.integers(6,size=13)]
  if len(set(order))!=6 or any(a==b for a,b in zip(order,order[1:])):continue
  schedule=tuple(e if i%2==0 else e[::-1] for i,e in enumerate(order))
  if key(schedule) in seen or not admissible(schedule):continue
  schedules.append(schedule);seen.add(key(schedule))
 assert len(schedules)==12
 plan.write_text(json.dumps({'seed':SEED,'schedules':schedules,'excluded_prior_files':coverage,'scope':'root search13 native candidate files only; not all previous experiments','short_maxiter':120,'deep_maxiter':500,'deep_count':3,'initialization':'normalized independent four-dimensional Gaussian quaternion for each SU2 local gate'},indent=2)+'\n')
 torch.set_num_threads(1);fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];v=torch.as_tensor(right_shift(4),dtype=torch.complex128);short=[];xs={}
 for case,schedule in enumerate(schedules):
  seed=SEED+100+case;x0=haar(np.random.default_rng(seed));x,record=fit(schedule,x0,fixed,v,120);xs[case]=x
  row={'case':case,'seed':seed,'schedule':schedule,**record,**save(case,'short',schedule,x)};short.append(row);print(json.dumps(row),flush=True)
 deep=[]
 if not any(r['valid_diagonalizer'] for r in short):
  for row in sorted(short,key=lambda r:(r['loss'],r['case']))[:3]:
   case=row['case'];x,record=fit(schedules[case],xs[case],fixed,v,500)
   r={'case':case,'schedule':schedules[case],'short_candidate':row['candidate'],**record,**save(case,'deep',schedules[case],x)};deep.append(r);print(json.dumps(r),flush=True)
 best=min(short+deep,key=lambda r:r['off_diagonal_error'])
 out.write_text(json.dumps({'plan':plan.name,'short_rows':short,'deep_rows':deep,'best_candidate':best['candidate'],'scope':'12 specified fresh schedules, one Haar local start each, three local continuations; no topology exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
