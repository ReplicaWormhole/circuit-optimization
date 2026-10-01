"""Owned differentiable XX/YY kernel and native-to-two-CX compiler."""
import math
import torch

def f_local(a,b):
    zero=a*0;even=a-b;odd=a+b
    ce,co=torch.cos(even),torch.cos(odd);se,so=1j*torch.sin(even),1j*torch.sin(odd)
    return torch.stack([torch.stack([ce,zero,zero,se]),torch.stack([zero,co,so,zero]),torch.stack([zero,so,co,zero]),torch.stack([se,zero,zero,ce])])

def f_matrix(a,b,first,second,n=4):
    assert first!=second and first in range(n) and second in range(n)
    local=f_local(a,b);indices=torch.arange(1<<n,device=a.device);mask1=1<<(n-1-first);mask2=1<<(n-1-second)
    bits=2*((indices&mask1)!=0).long()+((indices&mask2)!=0).long();rest=indices&~(mask1|mask2)
    return local[bits[:,None],bits[None,:]]*(rest[:,None]==rest[None,:])

def native_gate(first,second,a,b):
    return {'gate':'xx_yy','qubits':[first,second],'a':float(a),'b':float(b)}

def compile_f(g):
    c,t=g['qubits'];a,b=float(g['a']),float(g['b'])
    return [{'gate':'rx','qubit':c,'theta':-math.pi/2},{'gate':'rx','qubit':t,'theta':-math.pi/2},
            {'gate':'cx','control':c,'target':t},{'gate':'rx','qubit':c,'theta':-2*a},
            {'gate':'rz','qubit':t,'theta':-2*b},{'gate':'cx','control':c,'target':t},
            {'gate':'rx','qubit':c,'theta':math.pi/2},{'gate':'rx','qubit':t,'theta':math.pi/2}]

def compile_native(candidate):
    gates=[]
    for g in candidate['gates']:
        gates.extend(compile_f(g) if g['gate']=='xx_yy' else [dict(g)])
    return {'n':candidate['n'],'gates':gates}
