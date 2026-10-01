"""Push commuting local rotations forward; preserve counted entanglers.

No optimized or exact-certified candidate is asserted. Final output Rz
phases can be removed because they commute with the diagonal target.
"""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
import numpy as np
from check_circuit import one_qubit_matrix,evaluate
from search13_fulltail_invariant_gauge_rank import circuit_matrix
HERE=Path(__file__).resolve().parent
SOURCE=HERE/'labelswap_svd_refine_candidate.json'
def rot(axis,x):return one_qubit_matrix({'gate':'r'+axis,'qubit':0,'theta':float(x)})
def euler(m):
 a=m[0,0];b=-m[0,1];theta=2*np.arctan2(abs(b),abs(a))
 if abs(b)<1e-12:return -2*np.angle(a),0.,0.
 if abs(a)<1e-12:return -2*np.angle(b),np.pi,0.
 return -np.angle(a)-np.angle(b),theta,np.angle(b)-np.angle(a)
def main():
 out=HERE/'canonical13_result.json';assert not out.exists();c=json.loads(SOURCE.read_text());layers=[[np.eye(2,dtype=complex) for _ in range(4)]];schedule=[]
 for g in c['gates']:
  if g['gate']=='cx':schedule.append((g['control'],g['target']));layers.append([np.eye(2,dtype=complex) for _ in range(4)])
  else:
   m=one_qubit_matrix(g);m/=np.sqrt(np.linalg.det(m));wire=g['qubit'];layers[-1][wire]=m@layers[-1][wire]
 assert len(schedule)==13 and len(layers)==14
 gates=[];angles=[];b=rot('y',np.pi/2)
 for slot,(control,target) in enumerate(schedule):
  for wire in range(4):
   m=layers[slot][wire]
   if wire not in (control,target):passed=m
   else:
    axis='z' if wire==control else 'x';m0=m if axis=='z' else b.conj().T@m@b;alpha,beta,gamma=euler(m0);remain=rot('y',beta)@rot(axis,gamma);passed=rot(axis,alpha)
    assert np.max(abs(passed@remain-m))<1e-12
    gates.extend([{'gate':'r'+axis,'qubit':wire,'theta':float(gamma)},{'gate':'ry','qubit':wire,'theta':float(beta)}]);angles.extend([float(gamma),float(beta)])
   layers[slot+1][wire]=layers[slot+1][wire]@passed
  gates.append({'gate':'cx','control':control,'target':target})
 final=[]
 for wire,m in enumerate(layers[-1]):
  alpha,beta,gamma=euler(m);assert np.max(abs(rot('z',alpha)@rot('y',beta)@rot('z',gamma)-m))<1e-12
  final.append({'wire':wire,'alpha':float(alpha),'beta':float(beta),'gamma':float(gamma)})
  gates.extend([{'gate':'rz','qubit':wire,'theta':float(gamma)},{'gate':'ry','qubit':wire,'theta':float(beta)},{'gate':'rz','qubit':wire,'theta':float(alpha)}]);angles.extend([float(gamma),float(beta),float(alpha)])
 canonical={'n':4,'gates':gates};u=circuit_matrix(c['gates']);w=circuit_matrix(gates);phase=np.trace(u.conj().T@w)/16;error=float(np.max(abs(w-phase*u)));assert error<1e-12
 stripped={'n':4,'gates':[g for i,g in enumerate(gates) if not (i>=65 and (i-65)%3==2)]};assert len(stripped['gates'])==73
 for name,candidate in [('canonical13_candidate.json',canonical),('canonical13_stripped_candidate.json',stripped)]:
  q=evaluate(candidate);assert q['cnot_count']==13 and q['valid_diagonalizer'];(HERE/name).write_text(json.dumps(candidate,indent=2)+'\n')
 report={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'equivalence_up_to_global_phase_error':error,'schedule':schedule,'canonical_angles':angles,'final_layer':final,'pi24_distances':[float(abs(x/np.pi*24-round(x/np.pi*24))) for x in angles],'scope':'Numerical canonicalization preserving13CNOT; not exact angle recognition or synthesis bound'};out.write_text(json.dumps(report,indent=2)+'\n');print('equivalence error',error,'angles/pi',np.round(np.array(angles)/np.pi,6).tolist())
if __name__=='__main__':main()
