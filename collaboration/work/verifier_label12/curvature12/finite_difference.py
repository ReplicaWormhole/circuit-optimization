"""Reserved independent two-step added-layer curvature diagnostic."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE.parents[0]/'parity12'));from check_candidates import matrix,embed,PAULI
sys.path.insert(0,str(HERE.parents[0]));from verify import cycle

def unitary(z):
 local=np.asarray(z[:108]).reshape(9,4,3);f=np.asarray(z[108:]).reshape(4,2);u=np.eye(16,dtype=complex);pairs=[(0,2),(1,3),(0,1),(2,3)]
 for layer in range(9):
  op=np.array([[1]],complex)
  for axis in local[layer]:op=np.kron(op,expm(-.5j*sum(axis[j]*PAULI[p] for j,p in enumerate(['rx','ry','rz']))))
  u=op@u
  if layer in (0,1,2,5):
   c,t=[(0,1),(2,1),(3,1)][layer] if layer in (0,1,2) else (1,2);u=matrix({'gates':[{'gate':'cx','control':c,'target':t}]})@u
  elif layer in (3,4,6,7):
   i={3:0,4:1,6:2,7:3}[layer];a,b=f[i];u=embed(expm(1j*a*np.kron(PAULI['rx'],PAULI['rx'])+1j*b*np.kron(PAULI['ry'],PAULI['ry'])),*pairs[i])@u
 return u

def objective(z):
 u=unitary(z);b=u@cycle()@u.conj().T;off=b-np.diag(np.diag(b));return float(np.sum(abs(off)**2)/16)
def main():
 source=Path(sys.argv[1]);audit=json.loads(source.read_text());point=np.asarray(audit['point116']);direction24=np.asarray(audit['least_eigenvector24']);assert point.shape==(116,) and direction24.shape==(24,);assert abs(np.linalg.norm(direction24)-1)<1e-10
 direction=np.zeros(116);direction[12:36]=direction24;base=objective(point);checks=[]
 for h in (1e-3,5e-4):
  minus,plus=objective(point-h*direction),objective(point+h*direction);checks.append({'h':h,'minus_loss':minus,'base_loss':base,'plus_loss':plus,'gradient':(plus-minus)/(2*h),'curvature':(plus-2*base+minus)/(h*h)})
 result={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'direction_norm':float(np.linalg.norm(direction)),'mapping':'24CorderABaxis coordinates embedded116[12:36]','checks':checks,'threshold':-1e-6,'both_negative_below_threshold':all(row['curvature']<-1e-6 for row in checks),'scope':'Two predeclared reserved diagnosticsteps; no optimization or exactproof.'}
 Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
