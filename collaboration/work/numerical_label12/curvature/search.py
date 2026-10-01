"""Bounded label-free refinement and one AD Hessian, numerical evidence only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from check_circuit import evaluate,right_shift
from delete14_search import u3_to_axis,axis_to_u3
from search13_adaptive_topology import matrix_for_schedule
from search13_fresh_chain_dynamic_rows import unitary
HERE=Path(__file__).resolve().parent

def main():
    cfg=json.loads((HERE/'config.json').read_text());assert not (HERE/'result.json').exists()
    source=ROOT/cfg['source'];assert hashlib.sha256(source.read_bytes()).hexdigest()==cfg['source_sha256']
    c=json.loads(source.read_text());locals_=[g for g in c['gates'] if g['gate']=='u3']
    assert len(locals_)==56 and [g['qubit'] for g in locals_]==list(range(4))*14
    schedule=cfg['optimizer_schedule'];assert [(g['control'],g['target']) for g in c['gates'] if g['gate']=='cx']==[tuple(e) for e in schedule if e]
    x=np.array([u3_to_axis(g) for g in locals_]).ravel();torch.set_num_threads(1)
    cx=[torch.eye(16,dtype=torch.complex128) if e is None else matrix_for_schedule([e])[0] for e in schedule]
    v=torch.tensor(right_shift(4),dtype=torch.complex128);begun=time.monotonic();deadline=begun+120;rows=[];curvature=None
    def guard():
        if time.monotonic()>deadline:raise TimeoutError('overall compute cap')
    def loss(z):
        u=unitary(z.reshape(14,4,3),cx);b=u@v@u.conj().T;off=b-torch.diag(torch.diagonal(b))
        return (off.abs()**2).sum().real/16
    def fit(z,maxiter):
        last=z.copy()
        def f(w):
            nonlocal last
            guard();t=torch.tensor(w,dtype=torch.float64,requires_grad=True);value=loss(t);g=torch.autograd.grad(value,t)[0]
            last=w.copy();return float(value.detach()),g.detach().numpy()
        try:
            r=minimize(f,z,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'maxfun':20*maxiter,'maxls':20,'ftol':1e-15,'gtol':1e-11})
            return r.x,{'nit':int(r.nit),'nfev':int(r.nfev),'loss':float(r.fun),'message':str(r.message),'timedout':False}
        except TimeoutError:return last,{'timedout':True}
    def save(name,z,meta):
        gates=[]
        for slot in range(14):
            for q in range(4):gates.append({'gate':'u3','qubit':q,**axis_to_u3(*z.reshape(14,4,3)[slot,q])})
            if slot<13 and schedule[slot] is not None:
                a,b=schedule[slot];gates.append({'gate':'cx','control':a,'target':b})
        p=HERE/(name+'.json');candidate={'n':4,'gates':gates};p.write_text(json.dumps(candidate,indent=2)+'\n')
        check=evaluate(candidate);assert check['cnot_count']==12
        row={'name':name,'axis_checkpoint':z.tolist(),'candidate':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'optimizer':meta,'checker':check};rows.append(row)
        print(json.dumps({'name':name,'optimizer':meta,'checker':check}),flush=True)
    with threadpool_limits(limits=1):
        x,meta=fit(x,100);save('base',x,meta)
        if not meta['timedout']:
            try:
                guard();t=torch.tensor(x,dtype=torch.float64,requires_grad=True);value=loss(t);gradient=torch.autograd.grad(value,t)[0].detach().numpy()
                h=torch.autograd.functional.hessian(loss,t,vectorize=True).detach().numpy();guard()
                symmetry=float(np.max(abs(h-h.T)));hs=(h+h.T)/2;values,vectors=np.linalg.eigh(hs);direction=vectors[:,0]
                k=np.argmax(abs(direction))
                if direction[k]<0:direction=-direction
                residual=float(np.linalg.norm(hs@direction-values[0]*direction))
                curvature={'base_axis':x.tolist(),'gradient_norm':float(np.linalg.norm(gradient)),'gradient':gradient.tolist(),'hessian':h.tolist(),'symmetry_max_error':symmetry,'eigenvalues':values.tolist(),'least_eigenvector':direction.tolist(),'least_eigenpair_residual':residual,'negative_threshold':-1e-7}
                (HERE/'curvature.json').write_text(json.dumps(curvature,indent=2)+'\n')
                print(json.dumps({k:v for k,v in curvature.items() if k in ['gradient_norm','symmetry_max_error','least_eigenpair_residual','eigenvalues']}),flush=True)
                if values[0]<-1e-7:
                    for sign in (1,-1):
                        guard();z,info=fit(x+sign*.1*direction,150);save('plus' if sign==1 else 'minus',z,info)
                        if info['timedout']:break
            except TimeoutError:curvature={'timedout':True}
    (HERE/'result.json').write_text(json.dumps({'config':cfg,'elapsed_seconds':time.monotonic()-begun,'rows':rows,'curvature_file':'curvature.json' if (HERE/'curvature.json').exists() else None,'curvature_timedout':bool(curvature and curvature.get('timedout')),'scope':'Finite floating Hessian/local search, no topology exclusion or exact certificate.'},indent=2)+'\n')
if __name__=='__main__':main()
