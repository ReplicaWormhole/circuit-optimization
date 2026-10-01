"""One reserved, finite fresh8 full98 eigenlabel multistart batch. No exactness claims."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import importlib.util
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import scipy
import torch
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    rid = int(HERE.name)
    cfg = json.loads((HERE/'config.json').read_text())
    if (HERE/'result.json').exists():
        raise RuntimeError('Refusing to repeat a reserved batch')
    assert digest(__file__) == cfg['code_sha256']
    for p,h in cfg['dependencies'].items():
        assert digest(ROOT/p)==h, p
    label_spec = importlib.util.spec_from_file_location('label_objective', ROOT/cfg['assignment_module'])
    label_module = importlib.util.module_from_spec(label_spec)
    label_spec.loader.exec_module(label_module)
    assert scipy.__version__ == cfg['scipy_version']
    spec = importlib.util.spec_from_file_location('full98', ROOT/cfg['family_module'])
    family = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(family)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    start = time.monotonic()
    deadline = start + cfg['internal_seconds']
    rows = []
    best_index = None
    stop = 'all_starts_finished'

    def save():
        out = {'run_id':rid, 'config':cfg, 'rows':rows, 'best_index':best_index,
               'elapsed_seconds':time.monotonic()-start, 'stop_reason':stop,
               'versions':{'numpy':np.__version__, 'scipy':scipy.__version__, 'torch':torch.__version__},
               'scope':'Bounded numerical batch; not a family exclusion or exact certificate'}
        tmp = HERE/'result.tmp'
        tmp.write_text(json.dumps(out,indent=2)+'\n')
        tmp.replace(HERE/'result.json')

    with threadpool_limits(limits=1):
        for index,seed in enumerate(cfg['seeds']):
            if time.monotonic()>=deadline:
                stop='overall_deadline_before_start'; break
            rng=np.random.default_rng(seed)
            sigma=cfg['local_sigmas'][index%len(cfg['local_sigmas'])]
            initial=np.r_[rng.normal(0,sigma,90),rng.uniform(-np.pi/2,np.pi/2,8)]
            best={'loss':float('inf'),'x':None,'grad':None}
            calls=0
            fit_start=time.monotonic()
            termination='optimizer_returned'
            history=[]
            opt=None
            class BudgetStop(Exception): pass
            class NearHit(Exception): pass
            def fun(z):
                nonlocal calls,best
                if time.monotonic()>=deadline or calls>=cfg['maxfun']:
                    raise BudgetStop()
                x=torch.tensor(z,dtype=torch.float64,requires_grad=True)
                u=family.matrix(x)
                v_np=family.right_shift(4).astype(complex)
                labels,columns,cost=label_module.assignment(u.detach().numpy(),v_np)
                d=torch.tensor(labels,dtype=torch.complex128)
                residual=u@torch.tensor(v_np,dtype=torch.complex128)-d[:,None]*u
                value=(residual.abs()**2).sum().real/16
                grad=torch.autograd.grad(value,x)[0].detach().numpy()
                loss=float(value.detach())
                if not np.isfinite(loss) or not np.isfinite(grad).all():
                    raise FloatingPointError('Nonfinite fit evaluation')
                assert abs(loss-cost)<1e-12
                calls+=1
                history.append({'evaluation':calls,'loss':loss,'columns':columns.tolist()})
                if loss<best['loss']:
                    best={'loss':loss,'x':np.array(z).copy(),'grad':grad.copy()}
                if loss<cfg['nearhit_loss']:
                    raise NearHit()
                return loss,grad
            try:
                opt=minimize(fun,initial,jac=True,method='L-BFGS-B',options={
                    'maxiter':cfg['maxiter'],'maxfun':cfg['maxfun'],'maxls':cfg['maxls'],
                    'ftol':1e-16,'gtol':1e-12})
            except BudgetStop:
                termination='time_or_evaluation_budget'
            except NearHit:
                termination='promising_nearhit_requires_independent_verification'
            except FloatingPointError:
                termination='nonfinite_evaluation'
            if best['x'] is None:
                stop='no_evaluated_point'; break
            native=family.serialize(best['x'])
            compiled=family.helper.compile_native(native)
            npth=HERE/f'endpoint_{seed}_native.json'
            cpth=HERE/f'endpoint_{seed}_compiled.json'
            npth.write_text(json.dumps(native,indent=2)+'\n')
            cpth.write_text(json.dumps(compiled,indent=2)+'\n')
            check=family.evaluate(compiled)
            row={'seed':seed,'sigma':sigma,'initial':initial.tolist(), 'best':best['x'].tolist(),
                 'loss':best['loss'],'original_loss':float(family.objective(torch.tensor(best['x'],dtype=torch.float64)).detach()),'history':history,'gradient_norm':float(np.linalg.norm(best['grad'])),
                 'evaluations':calls,'elapsed_seconds':time.monotonic()-fit_start,
                 'termination':termination,'nit':None if opt is None else int(opt.nit),
                 'optimizer_message':None if opt is None else str(opt.message),
                 'native':npth.name,'compiled':cpth.name,
                 'native_sha256':digest(npth),'compiled_sha256':digest(cpth),'checker':check}
            rows.append(row)
            if best_index is None or row['loss']<rows[best_index]['loss']:
                best_index=len(rows)-1
                (HERE/'best_native.json').write_text(npth.read_text())
                (HERE/'best_compiled.json').write_text(cpth.read_text())
            print(json.dumps({k:row[k] for k in ('seed','sigma','loss','gradient_norm','evaluations','nit','termination')}),flush=True)
            save()
            if row['loss']<cfg['nearhit_loss']:
                stop='promising_nearhit_requires_independent_verification'; break
            if termination=='time_or_evaluation_budget' and time.monotonic()>=deadline:
                stop='overall_deadline'; break
    save()
    print(json.dumps({'run_id':rid,'completed_starts':len(rows),'best_loss':None if best_index is None else rows[best_index]['loss'], 'stop_reason':stop}),flush=True)

if __name__=='__main__':
    main()
