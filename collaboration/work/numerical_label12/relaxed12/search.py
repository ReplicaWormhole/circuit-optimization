"""Two embedded restricted-source relaxed-prefix fits; no exactness claim."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,hashlib,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from delete14_search import axis_rotation,axis_to_u3,kron_gate
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import matrix_for_schedule
from native_gate_check import evaluate as native_evaluate
HERE=Path(__file__).resolve().parent
helperpath=ROOT/'collaboration/work/verifier_label12/parity12/native_helpers.py'
spec=importlib.util.spec_from_file_location('parity12_helper',helperpath);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
def matrix(z):
    local=z[:108].reshape(9,4,3);f=z[108:].reshape(4,2);u=torch.eye(16,dtype=torch.complex128)
    pairs=[(0,2),(1,3),(0,1),(2,3)];cx=matrix_for_schedule([(0,1),(2,1),(3,1),(1,2)])
    for layer in range(9):
        for q in range(4):u=kron_gate(axis_rotation(*local[layer,q]),q)@u
        if layer in (0,1,2):u=cx[layer]@u
        elif layer==5:u=cx[3]@u
        elif layer in (3,4,6,7):
            i={3:0,4:1,6:2,7:3}[layer];a,b=pairs[i];u=helper.f_matrix(f[i,0],f[i,1],a,b)@u
    return u
def serialize(z):
    local=z[:108].reshape(9,4,3);f=z[108:].reshape(4,2);gates=[];pairs=[(0,2),(1,3),(0,1),(2,3)]
    for layer in range(9):
        for q in range(4):gates.append({'gate':'u3','qubit':q,**axis_to_u3(*local[layer,q])})
        if layer in (0,1,2):
            a,b=[(0,1),(2,1),(3,1)][layer];gates.append({'gate':'cx','control':a,'target':b})
        elif layer==5:gates.append({'gate':'cx','control':1,'target':2})
        elif layer in (3,4,6,7):
            i={3:0,4:1,6:2,7:3}[layer];a,b=pairs[i];gates.append(helper.native_gate(a,b,*f[i]))
    return {'n':4,'gates':gates}
def main():
    cfg=json.loads((HERE/'config.json').read_text());assert not (HERE/'result.json').exists()
    for name,digest in cfg['dependencies'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((ROOT/cfg['source']).read_bytes()).hexdigest()==cfg['source_sha256']
    torch.set_num_threads(1);v=torch.tensor(right_shift(4),dtype=torch.complex128);begun=time.monotonic();deadline=begun+120;rows=[]
    with threadpool_limits(limits=1):
        for entry in cfg['warmstarts']:
            seed=entry['source_seed'];initial=np.asarray(entry['params116'],dtype=float);assert initial.shape==(116,);x=initial.copy();meta={};timedout=False;initial_diagnostic={}
            def fun(z):
                nonlocal x
                if time.monotonic()>deadline:raise TimeoutError('overallcap')
                t=torch.tensor(z,dtype=torch.float64,requires_grad=True);u=matrix(t);b=u@v@u.conj().T;off=b-torch.diag(torch.diagonal(b));loss=(off.abs()**2).sum().real/16;g=torch.autograd.grad(loss,t)[0];x=z.copy()
                if not initial_diagnostic:
                    initial_diagnostic.update(loss=float(loss.detach()),gradient_norm=float(torch.linalg.vector_norm(g)),added_AB_gradient_norm=float(torch.linalg.vector_norm(g[12:36])))
                return float(loss.detach()),g.detach().numpy()
            try:
                r=minimize(fun,x,jac=True,method='L-BFGS-B',options={'maxiter':200,'maxfun':10000,'maxls':30,'ftol':1e-15,'gtol':1e-11});x=r.x;meta={'nit':int(r.nit),'nfev':int(r.nfev),'loss':float(r.fun),'message':str(r.message)}
            except TimeoutError:timedout=True
            native=serialize(x);compiled=helper.compile_native(native);p=HERE/f'native_seed{seed}.json';q=HERE/f'compiled_seed{seed}.json';p.write_text(json.dumps(native,indent=2)+'\n');q.write_text(json.dumps(compiled,indent=2)+'\n')
            ncheck=native_evaluate(native);ccheck=evaluate(compiled);assert ccheck['cnot_count']==12
            row={'seed':seed,'initial_diagnostic':initial_diagnostic,'initial_params':initial.tolist(),'params':x.tolist(),'optimizer':meta,'timedout':timedout,'native':str(p.relative_to(ROOT)),'compiled':str(q.relative_to(ROOT)),'native_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'compiled_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'native_checker':ncheck,'compiled_checker':ccheck};rows.append(row);print(json.dumps({k:w for k,w in row.items() if k not in ['initial_params','params']}),flush=True)
            if timedout:break
    (HERE/'result.json').write_text(json.dumps({'config':cfg,'elapsed_seconds':time.monotonic()-begun,'rows':rows,'scope':'Finite numerical family fits, passing candidate needs independent exact certification; no topology exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
