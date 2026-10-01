"""Reduced v=0 family: five free x/y SU2 layers, computational u preserved.

Import constructs no matrix/tensor/RNG. n3=(x,u,y), n4=(x,u,y,v), MSB first.
The reduced target is Wv0=CZ_xy CX(u->x), not a cyclic shift.
"""
from __future__ import annotations
import math
import numpy as np
import torch


def cx_numpy(control,target,n=3):
    out=np.zeros((1<<n,1<<n),dtype=np.complex128)
    for col in range(1<<n):
        row=col^(1<<(n-1-target)) if col&(1<<(n-1-control)) else col
        out[row,col]=1
    return out


def target_v0_numpy():
    cz=np.diag([-1 if (col&4) and (col&1) else 1 for col in range(8)]).astype(np.complex128)
    return cz@cx_numpy(1,0)


def rotvec_numpy(r):
    r=np.asarray(r,dtype=float);length=float(np.linalg.norm(r))
    factor=math.sin(length/2)/length if length else .5
    x,y,z=r*factor;c=math.cos(length/2)
    return np.array([[c-1j*z,-y-1j*x],[y-1j*x,c+1j*z]],dtype=np.complex128)


def matrix_rotvec(U):
    U=np.asarray(U,dtype=np.complex128);Q=U/np.sqrt(np.linalg.det(U))
    if float(np.trace(Q).real)<0:Q=-Q
    scalar=float(np.trace(Q).real/2)
    vector=np.array([-(Q[0,1]+Q[1,0]).imag/2,(Q[1,0]-Q[0,1]).real/2,-(Q[0,0]-Q[1,1]).imag/2])
    length=float(np.linalg.norm(vector))
    return np.zeros(3) if length<1e-14 else vector*(2*math.atan2(length,scalar)/length)


def collapse(source,checker):
    """Collapse the OLD source; source angle evaluation occurs only on call."""
    if source['n']!=3:raise ValueError('source must use n3=(x,u,y)')
    layers=[];current=[np.eye(2,dtype=np.complex128) for _ in range(2)];topology=[]
    for gate in source['gates']:
        if gate['gate']=='cx':
            layers.append(current);current=[np.eye(2,dtype=np.complex128) for _ in range(2)]
            topology.append([gate['control'],gate['target']])
        else:
            if gate['qubit'] not in (0,2):raise ValueError('all computational-u locals must be absent/identity; source contains u local')
            wire=0 if gate['qubit']==0 else 1
            current[wire]=checker.one_qubit_matrix(gate)@current[wire]
    layers.append(current)
    if topology!=[[0,2],[0,2],[1,0],[1,2]] or len(layers)!=5:raise ValueError('phase-correct source OLD order must be xy,xy,ux,uy')
    base=np.array([[matrix_rotvec(U) for U in layer] for layer in layers]).reshape(30)
    return base,topology,layers


def reduced_numpy(vectors,topology):
    layers=np.asarray(vectors).reshape(5,2,3);U=np.eye(8,dtype=np.complex128)
    for k,layer in enumerate(layers):
        M=np.kron(np.kron(rotvec_numpy(layer[0]),np.eye(2)),rotvec_numpy(layer[1]));U=M@U
        if k<4:U=cx_numpy(*topology[k])@U
    return U


def u3_from_su2(U,qubit):
    a,b,c,d=U[0,0],U[0,1],U[1,0],U[1,1];theta=2*math.atan2(abs(c),abs(a))
    if abs(a)>1e-12:
        phase=np.angle(a);phi=np.angle(c)-phase if abs(c)>1e-12 else 0.
        lam=np.angle(-b)-phase if abs(b)>1e-12 else np.angle(d)-phase
    else:phase=np.angle(c);phi=0.;lam=np.angle(-b)-phase
    return {'gate':'u3','qubit':qubit,'theta':float(theta),'phi':float(phi),'lam':float(lam)}


def serialize_reduced(vectors,topology):
    gates=[]
    for k,layer in enumerate(np.asarray(vectors).reshape(5,2,3)):
        gates.extend(u3_from_su2(rotvec_numpy(layer[j]),q) for j,q in enumerate((0,2)))
        if k<4:gates.append({'gate':'cx','control':topology[k][0],'target':topology[k][1]})
    return {'n':3,'gates':gates,'metadata':{'description':'five free x/y SU2 layers with computational u preserved','target':'Wv0=CZxy CXu-to-x; not V3','evidence_scope':'bounded reduced numerical diagnostic'}}


def serialize_full(vectors,topology,prefix):
    if prefix['n']!=4 or len(prefix['gates'])!=6:raise ValueError('accepted411 six-gate prefix required')
    gates=[dict(gate) for gate in prefix['gates']]
    for k,layer in enumerate(np.asarray(vectors).reshape(5,2,3)):
        gates.extend(u3_from_su2(rotvec_numpy(layer[j]),q) for j,q in enumerate((0,2)))
        if k<4:
            gates.append({'gate':'cx','control':topology[k][0],'target':topology[k][1]})
            if k==0:gates.append({'gate':'cx','control':3,'target':2})
    gates.append({'gate':'cx','control':3,'target':2})
    return {'n':4,'gates':gates,'metadata':{'description':'accepted411 Bell/router prefix plus reduced v0 new suffix and two bare CX32 comparators','evidence_scope':'full V4 numerical comparator; v1 extension unknown; reduced success is not full success'}}


class TorchWord:
    def __init__(self,topology):
        self.pauli=torch.tensor([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=torch.complex128)
        self.eye2=torch.eye(2,dtype=torch.complex128)
        self.cx=[torch.tensor(cx_numpy(*edge),dtype=torch.complex128) for edge in topology]
        self.target=torch.tensor(target_v0_numpy(),dtype=torch.complex128)

    def unitary(self,x):
        r=x.reshape(5,2,3);lengths=torch.linalg.vector_norm(r,dim=-1)
        factors=.5*torch.sinc(lengths/(2*math.pi));axes=torch.sum(r[...,None,None]*self.pauli,dim=-3)
        locals_=torch.cos(lengths/2)[...,None,None]*self.eye2-1j*factors[...,None,None]*axes
        U=torch.eye(8,dtype=torch.complex128)
        for k,layer in enumerate(locals_):
            M=torch.kron(torch.kron(layer[0].contiguous(),self.eye2),layer[1].contiguous());U=M@U
            if k<4:U=self.cx[k]@U
        return U

    def loss(self,x):
        U=self.unitary(x);A=U@self.target@U.conj().T;off=A-torch.diag(torch.diagonal(A))
        return torch.sum(torch.abs(off)**2).real/8
