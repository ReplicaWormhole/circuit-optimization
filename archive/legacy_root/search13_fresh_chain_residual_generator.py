"""Numerical nearest-sector correction, not a counted circuit synthesis.

Qubit 0 is the leftmost Pauli factor. W from polar(sum Q_l P_l)
intertwines numerical spectral projectors of D with nearest-root row masks.
No arbitrary four-qubit correction is counted as a free operation.
"""
import hashlib,itertools,json
from pathlib import Path
import numpy as np
from scipy.linalg import polar,logm
from check_circuit import right_shift
from search13_fulltail_invariant_gauge_rank import circuit_matrix
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_refine_trial1.json'
OUT=ROOT/'search13_fresh_chain_residual_generator_result.json'

def main():
 assert not OUT.exists()
 u=circuit_matrix(json.loads(SOURCE.read_text())['gates']);v=right_shift(4);d=u@v@u.conj().T
 roots=np.array([1,1j,-1,-1j]);labels=np.argmin(abs(np.diag(d)[:,None]-roots[None,:]),axis=1)
 assert np.bincount(labels,minlength=4).tolist()==[6,3,4,3]
 projectors=[sum(root**(-k)*np.linalg.matrix_power(d,k) for k in range(4))/4 for root in roots]
 s=sum(np.diag(labels==j)@p for j,p in enumerate(projectors));w,_=polar(s)
 target=np.diag(roots[labels]);intertwining=np.max(abs(w@d@w.conj().T-target));assert intertwining<1e-12
 h=1j*logm(w);assert np.max(abs(h-h.conj().T))<1e-12
 paulis={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
 rows=[];reconstructed=np.zeros((16,16),complex)
 for chars in itertools.product('IXYZ',repeat=4):
  p=np.array([[1]],complex)
  for c in chars:p=np.kron(p,paulis[c])
  a=np.trace(p@h)/16;assert abs(a.imag)<1e-12
  reconstructed+=a.real*p
  rows.append({'pauli':''.join(chars),'coefficient':float(a.real),'weight':sum(c!='I' for c in chars)})
 assert np.max(abs(h-reconstructed))<1e-12
 rows.sort(key=lambda r:abs(r['coefficient']),reverse=True)
 result={'source':SOURCE.name,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'labels':labels.tolist(),'intertwining_max_error':float(intertwining),'minimum_overlap_singular_value':float(np.linalg.svd(s,compute_uv=False).min()),'unitarity_error':float(np.max(abs(w.conj().T@w-np.eye(16)))),'generator_reconstruction_error':float(np.max(abs(h-reconstructed))),'generator_frobenius_norm':float(np.linalg.norm(h)),'pauli_coefficients':rows,'squared_coefficients_by_weight':{str(k):sum(r['coefficient']**2 for r in rows if r['weight']==k) for k in range(5)},'scope':'Numerical spectral-sector correction W only; W has not been synthesized or assigned CNOT cost. No upper/lower bound.'}
 OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['intertwining_max_error','minimum_overlap_singular_value','squared_coefficients_by_weight']}));print(json.dumps(rows[:12]))
if __name__=='__main__':main()
