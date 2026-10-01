"""Independent NumPy/SciPy validation of every complete gate list in run 389."""
import hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'experiments/runs/389'

def cx(c,t):
    M=np.zeros((16,16),complex)
    for j in range(16): M[j^(1<<(3-t)) if ((j>>(3-c))&1) else j,j]=1
    return M
PAULI={'rx':np.array([[0,1],[1,0]],complex),'ry':np.array([[0,-1j],[1j,0]],complex),'rz':np.diag([1,-1]).astype(complex)}
def embed(A,q):
    out=np.array([[1]],complex)
    for k in range(4): out=np.kron(out,A if k==q else np.eye(2))
    return out
def gate_matrix(g):
    kind=g['gate']
    if kind=='cx': return cx(g['control'],g['target'])
    if kind=='xx_yy':
        c,t=g['qubits']; X=PAULI['rx'];Y=PAULI['ry']
        return expm(1j*(float(g['a'])*embed(X,c)@embed(X,t)+float(g['b'])*embed(Y,c)@embed(Y,t)))
    if kind=='u3':
        th,ph,la=map(float,(g['theta'],g['phi'],g['lam']))
        A=np.array([[np.cos(th/2),-np.exp(1j*la)*np.sin(th/2)],
                    [np.exp(1j*ph)*np.sin(th/2),np.exp(1j*(ph+la))*np.cos(th/2)]],complex)
        return embed(A,int(g['qubit']))
    if kind in PAULI: return embed(expm(-0.5j*float(g['theta'])*PAULI[kind]),int(g['qubit']))
    raise ValueError(f'unsupported gate {g}')
def circuit(path):
    data=json.loads(path.read_text()); U=np.eye(16,dtype=complex)
    for g in data['gates']: U=gate_matrix(g)@U
    return data,U
def err_phase(A,B):
    p=np.vdot(B,A);p/=abs(p)
    return float(np.max(np.abs(A-p*B)))
def right_shift():
    V=np.zeros((16,16),complex)
    for j in range(16): V[((j&1)<<3)|(j>>1),j]=1
    return V
native_p=RUN/'base_native.json';compiled_p=RUN/'base_compiled.json'
nd,Un=circuit(native_p);cd,Uc=circuit(compiled_p);V=right_shift()
B=Un@V@Un.conj().T;off=B-np.diag(np.diag(B));fro=float(np.linalg.norm(off))
loss=fro*fro/16
result=json.loads((RUN/'result.json').read_text())
H=np.asarray(result['hessian'],float);ev=np.asarray(result['eigenvalues'],float);Q=np.asarray(result['eigenvectors_columns'],float)
hsym=(H+H.T)/2
report={
 'run_id':389,
 'native_sha256':hashlib.sha256(native_p.read_bytes()).hexdigest(),
 'compiled_sha256':hashlib.sha256(compiled_p.read_bytes()).hexdigest(),
 'native_total_gates':len(nd['gates']), 'native_cx_count':sum(g['gate']=='cx' for g in nd['gates']),
 'native_xx_yy_count':sum(g['gate']=='xx_yy' for g in nd['gates']), 'compiled_total_gates':len(cd['gates']),
 'compiled_cx_count':sum(g['gate']=='cx' for g in cd['gates']),
 'native_compiled_matrix_error_up_to_phase':err_phase(Un,Uc),
 'native_unitarity_max_error':float(np.max(np.abs(Un.conj().T@Un-np.eye(16)))),
 'target_offdiagonal_frobenius':fro,'target_offdiagonal_max_entry':float(np.max(np.abs(off))),
 'independent_normalized_loss':loss,'saved_loss':float(result['loss']),
 'loss_abs_difference':abs(loss-float(result['loss'])),
 'independent_hessian_symmetry_max':float(np.max(np.abs(H-H.T))),
 'recomputed_least_eigenvalue':float(np.linalg.eigvalsh(hsym)[0]),
 'saved_least_eigenvalue':float(ev[0]),
 'least_eigenpair_residual':float(np.linalg.norm(hsym@Q[:,0]-ev[0]*Q[:,0])),
 'exact_certificate':False,
 'assessment':'Both saved gate lists are complete, unitary, and equivalent; the 12-CX list does not diagonalize V4. The saved AD Hessian is symmetric and its eigensystem is internally consistent; these checks do not independently verify derivatives.'}
assert report['compiled_cx_count']==12 and report['native_cx_count']==4 and report['native_xx_yy_count']==4
assert report['native_compiled_matrix_error_up_to_phase']<2e-12 and report['native_unitarity_max_error']<2e-12
assert report['loss_abs_difference']<1e-12 and report['least_eigenpair_residual']<1e-9
out=Path(__file__).with_name('run389_gate_audit.json');out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
