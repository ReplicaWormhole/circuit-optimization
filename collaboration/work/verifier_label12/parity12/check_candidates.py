"""Independent SciPy exponential/native and NumPy compiled matrix comparison."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[0]));from verify import matrix as u3_matrix,cycle
PAULI={'rx':np.array([[0,1],[1,0]]),'ry':np.array([[0,-1j],[1j,0]]),'rz':np.diag([1,-1])}
def embed(local,c,t):
 op=np.zeros((16,16),complex);m1,m2=1<<(3-c),1<<(3-t)
 for x in range(16):
  col=2*bool(x&m1)+bool(x&m2);rest=x&~(m1|m2)
  for row in range(4):op[rest|(m1 if row&2 else 0)|(m2 if row&1 else 0),x]=local[row,col]
 return op
def matrix(data):
 u=np.eye(16,dtype=complex)
 for g in data['gates']:
  if g['gate']=='xx_yy':
   c,t=g['qubits'];x,y=PAULI['rx'],PAULI['ry'];op=embed(expm(1j*g['a']*np.kron(x,x)+1j*g['b']*np.kron(y,y)),c,t)
  elif g['gate'] in PAULI:
   local=expm(-1j*g['theta']*PAULI[g['gate']]/2);op=np.array([[1]],complex)
   for q in range(4):op=np.kron(op,local if q==g['qubit'] else np.eye(2))
  else:op,_=u3_matrix({'n':4,'gates':[g]})
  u=op@u
 return u
if __name__=='__main__':
 rows=[]
 for nativepath,compiledpath in zip(sys.argv[2::2],sys.argv[3::2]):
  n,c=Path(nativepath),Path(compiledpath);nd,cd=json.loads(n.read_text()),json.loads(c.read_text());u,w=matrix(nd),matrix(cd);b=w@cycle()@w.conj().T;off=b-np.diag(np.diag(b));labels=np.argmin(abs(np.diag(b)[:,None]-np.array([1,1j,-1,-1j])),axis=1)
  row={'native':str(n),'compiled':str(c),'native_sha256':hashlib.sha256(n.read_bytes()).hexdigest(),'compiled_sha256':hashlib.sha256(c.read_bytes()).hexdigest(),'native_cx':sum(g['gate']=='cx' for g in nd['gates']),'native_f':sum(g['gate']=='xx_yy' for g in nd['gates']),'compiled_cx':sum(g['gate']=='cx' for g in cd['gates']),'full_matrix_difference':float(np.max(abs(u-w))),'maxoff':float(np.max(abs(off))),'offloss':float(np.sum(abs(off)**2)/16),'unitarity':float(np.max(abs(w.conj().T@w-np.eye(16)))),'nearest_labels':labels.tolist(),'root_counts':np.bincount(labels,minlength=4).tolist(),'exact_certificate':False}
  assert row['native_cx']==4 and row['native_f']==4 and row['compiled_cx']==12;assert row['full_matrix_difference']<1e-12;rows.append(row)
 Path(sys.argv[1]).write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
