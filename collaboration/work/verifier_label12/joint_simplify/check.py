"""Independent local-rule replay and NumPy matrix equivalence."""
import json,sys,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
folder=ROOT/'collaboration/work/numerical_label12/joint_simplify'
r=json.loads((folder/'result.json').read_text());source=ROOT/r['config']['source'];candidate=ROOT/r['candidate'];gates=json.loads(source.read_text())['gates']
def support(g):return {g['control'],g['target']} if g['gate']=='cx' else {g['qubit']}
def commute(a,b):
 if support(a).isdisjoint(support(b)):return True
 if a['gate']==b['gate']=='cx':return a['control']!=b['target'] and b['control']!=a['target']
 if a['gate']=='cx':return b['gate']=='rz' and b['qubit']==a['control'] or b['gate']=='rx' and b['qubit']==a['target']
 if b['gate']=='cx':return commute(b,a)
 return a['gate']==b['gate'] and a['gate'] in ['rx','ry','rz']
def rational(g):return Fraction(g['theta'].split('*pi')[0])/int(g['theta'].split('/')[-1])
for event in r['trace']:
 i,j=event['i'],event['j'];a,b=gates[i],gates[j];assert a==event['a'] and b==event['b'];assert gates[i+1:j]==event['crossed'];assert all(commute(a,m) for m in gates[i+1:j])
 if event['rule']=='involution cancellation':assert a==b and a['gate'] in ['h','cx'];assert not event['replacement']
 else:
  assert a['gate']==b['gate'] and a['qubit']==b['qubit'];total=rational(a)+rational(b);assert (not event['replacement'] and total==0) or rational(event['replacement'][0])==total
 gates=gates[:i]+event['replacement']+gates[i+1:j]+gates[j+1:]
assert gates==json.loads(candidate.read_text())['gates']
def matrix(gates):
 u=np.eye(16,dtype=complex)
 for g in gates:
  if g['gate']=='cx':
   op=np.zeros((16,16),complex)
   for x in range(16):op[x^(1<<(3-g['target'])) if (x>>(3-g['control']))&1 else x,x]=1
  else:
   if g['gate']=='h':local=np.array([[1,1],[1,-1]])/np.sqrt(2)
   elif g['gate']=='x':local=np.array([[0,1],[1,0]])
   else:
    t=float(rational(g))*np.pi;c,s=np.cos(t/2),np.sin(t/2)
    if g['gate']=='ry':local=np.array([[c,-s],[s,c]])
    elif g['gate']=='rz':local=np.diag([np.exp(-1j*t/2),np.exp(1j*t/2)])
    else:raise ValueError(g)
   op=np.array([[1]],complex)
   for q in range(4):op=np.kron(op,local if q==g['qubit'] else np.eye(2))
  u=op@u
 return u
u=matrix(json.loads(source.read_text())['gates']);w=matrix(gates);v=np.zeros((16,16),complex)
for x in range(16):v[(x>>1)|((x&1)<<3),x]=1
b=w@v@w.conj().T;off=b-np.diag(np.diag(b));report={'trace_exact_replay_verified':True,'rewrites':len(r['trace']),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'candidate_sha256':hashlib.sha256(candidate.read_bytes()).hexdigest(),'cx_count':sum(g['gate']=='cx' for g in gates),'maximum_full_matrix_difference':float(np.max(abs(u-w))),'maxoff':float(np.max(abs(off))),'unitarity_error':float(np.max(abs(w.conj().T@w-np.eye(16)))),'whole_exact_UV_certificate':False}
assert report['maximum_full_matrix_difference']<1e-12;assert report['maxoff']<1e-12
Path(__file__).with_name('independent_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
