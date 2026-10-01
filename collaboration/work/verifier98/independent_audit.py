"""Independent NumPy/SciPy audit of the 116-to-98 map and CX+F compilation.
No fitting or search is performed. Qubit 0 is the most-significant bit.
"""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parents[3]
P=ROOT/'collaboration/work/numerical_label12/relaxed12/config.json'
source=json.loads(P.read_text())
entry=next(e for e in source['warmstarts'] if e['source_seed']==9261101)
z=np.asarray(entry['params116'],dtype=float)
assert z.shape==(116,)

def su2(v):
    x,y,z=v; t=float(np.linalg.norm([x,y,z])); I=np.eye(2,dtype=complex)
    if t==0: return I
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
    return np.cos(t/2)*I-1j*np.sin(t/2)/t*(x*X+y*Y+z*Z)

def kron_local(vs):
    out=np.array([[1]],complex)
    for v in vs: out=np.kron(out,su2(v))
    return out

def cx(c,t):
    M=np.zeros((16,16),complex)
    for j in range(16):
        k=j^(1<<(3-t)) if ((j>>(3-c))&1) else j
        M[k,j]=1
    return M

def fmat(a,b,c,t):
    paulis=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
    def embed(p,q):
        out=np.array([[1]],complex)
        for k in range(4): out=np.kron(out,p if k==q else np.eye(2))
        return out
    return expm(1j*(a*embed(paulis[0],c)@embed(paulis[0],t)+b*embed(paulis[1],c)@embed(paulis[1],t)))

def full(z116):
    loc=z116[:108].reshape(9,4,3); fs=z116[108:].reshape(4,2)
    pairs=[(0,2),(1,3),(0,1),(2,3)]
    U=np.eye(16,dtype=complex)
    for layer in range(9):
        U=kron_local(loc[layer])@U
        if layer in (0,1,2): U=cx(*[(0,1),(2,1),(3,1)][layer])@U
        elif layer==5: U=cx(1,2)@U
        elif layer in (3,4,6,7):
            idx={3:0,4:1,6:2,7:3}[layer]; U=fmat(*fs[idx],*pairs[idx])@U
    return U

def vec_from_su2(M):
    # Remove roundoff-scale determinant phase, then recover principal SU(2) axis-angle.
    M=np.asarray(M,complex); M=M/np.sqrt(np.linalg.det(M))
    c=float(np.clip(np.trace(M).real/2,-1,1)); theta=2*math.acos(c)
    if theta<1e-12: return np.zeros(3)
    s=math.sin(theta/2)
    if abs(s)<1e-10: raise ValueError('singular SU2 chart at angle pi')
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
    return theta*np.array([np.trace(1j*M@X).real, np.trace(1j*M@Y).real, np.trace(1j*M@Z).real]) /(2*s)

def transform(z116):
    local=z116[:108].reshape(9,4,3); fs=z116[108:].reshape(4,2)
    mats=np.empty((9,4,2,2),complex)
    for l in range(9):
        for q in range(4): mats[l,q]=su2(local[l,q])
    # Original groups are L0,A,B,L1,L2,...,L6.
    new=[]
    L0,A,B,L1= mats[:4]
    n0=L0.copy(); n1=L1.copy()
    n0[2]=A[2]@L0[2]
    n0[3]=B[3]@A[3]@L0[3]
    n1[0]=L1[0]@B[0]@A[0]
    n1[2]=L1[2]@B[2]
    # Optimizer order: seven free layers [L0',L1',L2..L6], then A1,B1, then F.
    new.extend([[vec_from_su2(M) for M in n0], [vec_from_su2(M) for M in n1]])
    for l in range(4,9): new.append(local[l].copy())
    return np.concatenate([np.asarray(new).reshape(-1), local[1,1], local[2,1], fs.reshape(-1)])

def main():
    # A deterministic arbitrary nonzero point exercises every absorbed boundary product.
    rng=np.random.default_rng(140001); test=z.copy(); test[:108]=rng.normal(0,.31,108)
    zt=transform(test); assert zt.shape==(98,)
    # Rebuild 98 family in its literal chronology, directly from transformed SU2s.
    def full98(z98):
        free=z98[:84].reshape(7,4,3); aa=np.zeros((4,3));bb=np.zeros((4,3));aa[1]=z98[84:87];bb[1]=z98[87:90]
        loc=np.concatenate([free[:1],aa[None,:,:],bb[None,:,:],free[1:]],axis=0); fs=z98[90:].reshape(4,2)
        pairs=[(0,2),(1,3),(0,1),(2,3)];U=np.eye(16,dtype=complex)
        for layer in range(9):
            U=kron_local(loc[layer])@U
            if layer in (0,1,2): U=cx(*[(0,1),(2,1),(3,1)][layer])@U
            elif layer==5: U=cx(1,2)@U
            elif layer in (3,4,6,7):
                idx={3:0,4:1,6:2,7:3}[layer]; U=fmat(*fs[idx],*pairs[idx])@U
        return U
    U116=full(test); U98=full98(zt)
    phase=np.vdot(U116,U98);phase/=abs(phase)
    map_err=float(np.max(np.abs(U116-phase*U98)))
    # Zero-A/B exact embedding of frozen 385 seed: preserve original free boundaries.
    zero=z.copy(); zero[12:36]=0
    # explicit source mapping is identity on free layers/F at A=B=I
    # Construct 98 ordering: L0, A1=I, B1=I, L1,L2..L6,F.
    local=zero[:108].reshape(9,4,3)
    loc98=np.concatenate([local[0].ravel(),local[3].ravel(),* [local[k].ravel() for k in range(4,9)]])
    z98=np.concatenate([loc98,np.zeros(6),zero[108:]])
    embed_err=float(np.max(np.abs(full(zero)-full98(z98))))
    # Independent gate-list simulator for saved native/compiled outputs.
    def oneq(g):
        name=g['gate']
        if name=='u3':
            th,ph,la=g['theta'],g['phi'],g['lam'];return np.array([[np.cos(th/2),-np.exp(1j*la)*np.sin(th/2)],[np.exp(1j*ph)*np.sin(th/2),np.exp(1j*(ph+la))*np.cos(th/2)]],complex)
        t=g['theta']; P={'rx':np.array([[0,1],[1,0]],complex),'ry':np.array([[0,-1j],[1j,0]],complex),'rz':np.diag([1,-1]).astype(complex)}[name]
        return expm(-.5j*t*P)
    def saved_matrix(payload):
        data = payload if isinstance(payload, dict) else json.loads(Path(payload).read_text())
        gates=data['gates'];U=np.eye(16,dtype=complex)
        for g in gates:
            if g['gate']=='cx': G=cx(g['control'],g['target'])
            elif g['gate']=='xx_yy': G=fmat(g['a'],g['b'],*g['qubits'])
            elif g['gate'] in ('u3','rx','ry','rz'):
                q=g['qubit'];G=np.array([[1]],complex)
                for k in range(4):G=np.kron(G,oneq(g) if k==q else np.eye(2))
            else: raise ValueError(g)
            U=G@U
        return U,gates
    folder=ROOT/'collaboration/work/numerical_label12/relaxed12'
    native=folder/'native_seed9261101.json'; compiled=folder/'compiled_seed9261101.json'
    Un,gn=saved_matrix(native);Uc,gc=saved_matrix(compiled)
    ph=np.vdot(Un,Uc);ph/=abs(ph);comp_err=float(np.max(np.abs(Un-ph*Uc)))
    V=np.zeros((16,16),complex)
    for j in range(16): V[((j&1)<<3)|(j>>1),j]=1
    D=Un@V@Un.conj().T; off=D-np.diag(np.diag(D))
    report={
     'seed':9261101,'source_config_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),
     'params116_sha256':hashlib.sha256(np.asarray(z,dtype='<f8').tobytes()).hexdigest(),
     'params116_shape':list(z.shape),'full116_groups':'9 x 4 x 3 local-axis coordinates, then 4 x 2 F angles',
     '98_layout_verified':'optimizer: seven free layers [L0prime,L1prime,L2..L6] (84), A1q1(3), B1q1(3), F(8); chronology inserts A1,B1 between CX01/CX21/CX31',
     'random_nonzero_boundary_map_max_difference_up_to_phase':map_err,
     'frozen_seed_zero_AB_embedding_max_difference':embed_err,
     'native_gate_count':len(gn),'native_F_count':sum(g['gate']=='xx_yy' for g in gn),'native_CX_count':sum(g['gate']=='cx' for g in gn),
     'compiled_gate_count':len(gc),'compiled_CX_count':sum(g['gate']=='cx' for g in gc),
     'independent_native_to_compiled_max_difference_up_to_phase':comp_err,
     'compiled_cnot_count':sum(g['gate']=='cx' for g in gc),
     'independent_native_unitarity_error':float(np.max(np.abs(Un.conj().T@Un-np.eye(16)))),
     'independent_native_target_offdiagonal_frobenius':float(np.linalg.norm(off)),
     'native_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),'compiled_sha256':hashlib.sha256(compiled.read_bytes()).hexdigest(),
     'scope':'Structural boundary-map check, one frozen embedding and one saved gate-list audit; no optimization and no exact certificate.'}
    assert map_err<2e-12 and embed_err<2e-12 and comp_err<2e-12 and report['compiled_cnot_count']==12
    out=Path(__file__).with_name('independent_audit.json');out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__ == "__main__":
    main()
