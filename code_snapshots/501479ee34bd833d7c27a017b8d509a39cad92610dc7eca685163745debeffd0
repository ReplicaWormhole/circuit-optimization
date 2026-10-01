"""Two frozen unordered-pair mutations; numerical evidence only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
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
    p=ROOT/cfg['source'];assert hashlib.sha256(p.read_bytes()).hexdigest()==cfg['source_sha256'];c=json.loads(p.read_text())
    local=[g for g in c['gates'] if g['gate']=='u3'];assert len(local)==56 and [g['qubit'] for g in local]==list(range(4))*14
    start=np.array([u3_to_axis(g) for g in local]).ravel();torch.set_num_threads(1)
    v=torch.tensor(right_shift(4),dtype=torch.complex128);begun=time.monotonic();deadline=begun+120;rows=[]
    with threadpool_limits(limits=1):
        for index,proposal in enumerate(cfg['proposals']):
            schedule=proposal['schedule'];cx=[torch.eye(16,dtype=torch.complex128) if e is None else matrix_for_schedule([e])[0] for e in schedule];x=start.copy();timedout=False;meta={}
            def fun(z):
                nonlocal x
                if time.monotonic()>deadline:raise TimeoutError('overallcompute cap')
                t=torch.tensor(z,dtype=torch.float64,requires_grad=True);u=unitary(t.reshape(14,4,3),cx);b=u@v@u.conj().T;off=b-torch.diag(torch.diagonal(b));loss=(off.abs()**2).sum().real/16;gradient=torch.autograd.grad(loss,t)[0];x=z.copy()
                return float(loss.detach()),gradient.detach().numpy()
            try:
                r=minimize(fun,x,jac=True,method='L-BFGS-B',options={'maxiter':200,'maxfun':4000,'maxls':20,'ftol':1e-15,'gtol':1e-11});x=r.x;meta={'nit':int(r.nit),'nfev':int(r.nfev),'loss':float(r.fun),'message':str(r.message)}
            except TimeoutError:timedout=True
            gates=[]
            for slot in range(14):
                for q in range(4):gates.append({'gate':'u3','qubit':q,**axis_to_u3(*x.reshape(14,4,3)[slot,q])})
                if slot<13 and schedule[slot] is not None:
                    a,b=schedule[slot];gates.append({'gate':'cx','control':a,'target':b})
            candidate={'n':4,'gates':gates};p=HERE/f'target{index}.json';p.write_text(json.dumps(candidate,indent=2)+'\n');check=evaluate(candidate);assert check['cnot_count']==12
            row={'proposal':proposal,'candidate':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'axis_checkpoint':x.tolist(),'optimizer':meta,'timedout':timedout,'checker':check};rows.append(row);print(json.dumps({k:z for k,z in row.items() if k!='axis_checkpoint'}),flush=True)
            if timedout:break
    (HERE/'result.json').write_text(json.dumps({'config':cfg,'elapsed_seconds':time.monotonic()-begun,'rows':rows,'scope':'Finite two architecture fits, no topology exclusion or exact certificate'},indent=2)+'\n')
if __name__=='__main__':main()
