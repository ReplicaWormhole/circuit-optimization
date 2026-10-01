"""Frozen one-target spectral-label experiment; numerical evidence only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy.optimize import minimize,linear_sum_assignment
from threadpoolctl import threadpool_limits
from check_circuit import evaluate,right_shift
from delete14_search import u3_to_axis,axis_to_u3
from search13_adaptive_topology import matrix_for_schedule
from search13_fresh_chain_dynamic_rows import unitary
HERE=Path(__file__).resolve().parent

def main():
    cfg=json.loads((HERE/'config.json').read_text())
    assert not (HERE/'result.json').exists()
    source=ROOT/cfg['source'];assert hashlib.sha256(source.read_bytes()).hexdigest()==cfg['source_sha256']
    c=json.loads(source.read_text());schedule=cfg['optimizer_schedule']
    locals_=[g for g in c['gates'] if g['gate']=='u3']
    assert len(locals_)==56 and [g['qubit'] for g in locals_]==list(range(4))*14
    assert [(g['control'],g['target']) for g in c['gates'] if g['gate']=='cx']==[tuple(e) for e in schedule if e is not None]
    start=np.array([u3_to_axis(g) for g in locals_]).ravel()
    cx=[torch.eye(16,dtype=torch.complex128) if e is None else matrix_for_schedule([e])[0] for e in schedule]
    v=torch.tensor(right_shift(4),dtype=torch.complex128);roots=np.array([1,1j,-1,-1j])
    torch.set_num_threads(1)
    begun=time.monotonic();deadline=begun+cfg['wall_seconds'];rows=[]
    def guard():
        if time.monotonic()>deadline:raise TimeoutError('overall computation cap')
    def residual(z,labels):
        u=unitary(z.reshape(14,4,3),cx);d=torch.tensor(roots[labels],dtype=torch.complex128)
        r=(u@v-d[:,None]*u)/4
        return torch.cat((r.real.ravel(),r.imag.ravel()))
    def arr(z,labels):
        guard();return residual(torch.tensor(z,dtype=torch.float64),labels).detach().numpy()
    with threadpool_limits(limits=1):
        u=unitary(torch.tensor(start.reshape(14,4,3)),cx).detach().numpy();b=u@v.numpy()@u.conj().T
        repeated=np.repeat(np.arange(4),[6,3,4,3]);_,cols=linear_sum_assignment(abs(np.diag(b)[:,None]-roots[repeated][None,:])**2)
        labels=repeated[cols]
        r,s=10,14
        assert labels[r]!=labels[s]
        d=labels.copy();d[r],d[s]=d[s],d[r]
        targets=[{'pair':[r,s],'labels':d.tolist()}]
        (HERE/'targets.json').write_text(json.dumps({'base_labels':labels.tolist(),'targets':targets},indent=2)+'\n')
        for index,target in enumerate(targets):
            x=start.copy();d=np.array(target['labels']);trace=[];timedout=False;meta={}
            def fun(z):
                nonlocal x
                guard();t=torch.tensor(z,dtype=torch.float64,requires_grad=True);r=residual(t,d);loss=r@r
                grad=torch.autograd.grad(loss,t)[0];x=np.asarray(z).copy()
                return float(loss.detach()),grad.detach().numpy()
            try:
                opt=minimize(fun,x,jac=True,method='L-BFGS-B',options={'maxiter':200,'maxfun':4000,'maxls':20,'ftol':1e-15,'gtol':1e-11})
                x=opt.x.copy();meta={'nit':int(opt.nit),'nfev':int(opt.nfev),'message':str(opt.message),'loss':float(opt.fun)}
                for step_index in range(5):
                    guard();r=arr(x,d);old=float(r@r)
                    t=torch.tensor(x,dtype=torch.float64,requires_grad=True)
                    j=torch.autograd.functional.jacobian(lambda z:residual(z,d),t,vectorize=True).detach().numpy()
                    delta,_,rank,_=np.linalg.lstsq(j,-r,rcond=1e-6);accepted=False;new=old
                    for k in range(8):
                        trial=x+2.**(-k)*delta;rr=arr(trial,d);new=float(rr@rr)
                        if new<old:x=trial;accepted=True;break
                    trace.append({'step':step_index,'before':old,'after':new,'rank':int(rank),'accepted':accepted,'alpha':2.**(-k)})
                    if not accepted:break
            except TimeoutError:
                timedout=True
            gates=[]
            for slot in range(14):
                for q in range(4):gates.append({'gate':'u3','qubit':q,**axis_to_u3(*x.reshape(14,4,3)[slot,q])})
                if slot<13 and schedule[slot] is not None:
                    a,z=schedule[slot];gates.append({'gate':'cx','control':a,'target':z})
            candidate={'n':4,'gates':gates};path=HERE/f'target{index}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
            check=evaluate(candidate);assert check['cnot_count']==12
            row={**target,'optimizer':meta,'newton':trace,'timedout':timedout,'candidate':str(path.relative_to(ROOT)),'candidate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'checker':check}
            rows.append(row);print(json.dumps(row),flush=True)
            if timedout:break
    result={'config':cfg,'elapsed_seconds':time.monotonic()-begun,'rows':rows,'scope':'Finite numerical search; failure excludes no target family or topology; passing output needs exact certification.'}
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
