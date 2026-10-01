"""One frozen two-start, twelve-CX round; no exactness or exclusion claim."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import sys, json, time, hashlib
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from check_circuit import evaluate, right_shift
from delete14_search import axis_to_u3
from search13_adaptive_topology import matrix_for_schedule, objective
from search13_fresh_chain_dynamic_rows import decode
OUT = Path(__file__).resolve().parent
SOURCE = ROOT/'collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json'
ORIGINAL = [(0,1),(3,2),(1,2),(1,0),(2,3),(2,1),(0,1),(3,2),(1,2),(1,0),(2,3),(2,1),(0,1)]
SCHEDULES=[]
for motif in (0,6):
    schedule=ORIGINAL.copy();schedule[motif:motif+3]=[None,(1,2),(0,2)];SCHEDULES.append(schedule)
torch.set_num_threads(1)

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2)+'\n')

def matrix_u3(m, q):
    alpha = np.angle(m[0,0])
    theta = 2*np.arctan2(abs(m[1,0]), abs(m[0,0]))
    if abs(m[1,0]) < 1e-14:
        phi, lam = 0., np.angle(m[1,1])-alpha
    else:
        phi, lam = np.angle(m[1,0])-alpha, np.angle(-m[0,1])-alpha
    return dict(gate='u3', qubit=q, theta=float(theta), phi=float(phi), lam=float(lam))

def warm_source(source):
    paulis = {'rx':np.array([[0,1],[1,0]]), 'ry':np.array([[0,-1j],[1j,0]]), 'rz':np.diag([1,-1])}
    mats = [np.eye(2,dtype=complex) for _ in range(4)]
    gates = []
    for g in source['gates']:
        if g['gate']=='cx':
            gates.extend(matrix_u3(m,q) for q,m in enumerate(mats))
            gates.append(g.copy())
            mats = [np.eye(2,dtype=complex) for _ in range(4)]
        else:
            t=g['theta']; r=np.cos(t/2)*np.eye(2)-1j*np.sin(t/2)*paulis[g['gate']]
            mats[g['qubit']]=r@mats[g['qubit']]
    gates.extend(matrix_u3(m,q) for q,m in enumerate(mats))
    return {'n':4,'gates':gates}

def serialize(x, schedule):
    gates=[]
    for slot, layer in enumerate(np.asarray(x).reshape(14,4,3)):
        gates.extend(dict(gate='u3',qubit=q,**axis_to_u3(*a)) for q,a in enumerate(layer))
        if slot<13 and schedule[slot] is not None:
            c,t=schedule[slot];gates.append(dict(gate='cx',control=c,target=t))
    return {'n':4,'gates':gates,'evidence':'numerical diagnostic; requires independent exact certification'}

def independent(candidate):
    u=np.eye(16,dtype=complex)
    for g in candidate['gates']:
        if g['gate']=='cx':
            big=np.zeros((16,16),complex)
            for b in range(16):
                target=b^(1<<(3-g['target'])) if b&(1<<(3-g['control'])) else b
                big[target,b]=1
        else:
            t,p,l=g['theta'],g['phi'],g['lam'];c,s=np.cos(t/2),np.sin(t/2)
            small=np.array([[c,-np.exp(1j*l)*s],[np.exp(1j*p)*s,np.exp(1j*(p+l))*c]])
            big=np.array([[1]],complex)
            for q in range(4):big=np.kron(big,small if q==g['qubit'] else np.eye(2))
        u=big@u
    v=np.zeros((16,16),complex)
    for b in range(16):v[(b>>1)|((b&1)<<3),b]=1
    d=u@v@u.conj().T
    return {'off_diagonal_error':float(np.max(abs(d-np.diag(np.diag(d))))),'unitarity_error':float(np.max(abs(u@u.conj().T-np.eye(16)))), 'cnot_count':sum(g['gate']=='cx' for g in candidate['gates'])}

class BudgetExceeded(Exception):pass

def prepare():
    source=json.loads(SOURCE.read_text());temp=warm_source(source)
    assert evaluate(temp)['valid_diagonalizer']
    x=decode(temp).ravel()
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)]
    v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    original=[(g['control'],g['target']) for g in source['gates'] if g['gate']=='cx']
    assert objective(x,fixed,matrix_for_schedule(original),v,False)<1e-20
    for schedule in SCHEDULES:
        native=[s for s in schedule if s is not None]
        assert len(native)==12
        assert all(native != original[:k]+original[k+1:] for k in range(13))
    plan={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'schedules':SCHEDULES,'motif_slots':[0,6],'rule':'replace source01,32,12 triple by identity,12,02 motivated by decorated01,12,01 braid; refit all locals; no source equivalence asserted','seeds':[926301,926302],'perturbation_stds':[.08,.08],'optimizer':'SciPy L-BFGS-B analytic Torch gradient','maxiter_per_start':300,'per_start_wall_cap_seconds':120,'total_compute_wall_cap_seconds':240,'threads':1,'local_parameter_count':168,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependency_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['delete14_search.py','search13_adaptive_topology.py','search13_fresh_chain_dynamic_rows.py','check_circuit.py','collaboration/work/numerical_round2/run.py']}}
    return plan,x,fixed,v

def main():
    plan,x,fixed,v=prepare()
    if '--prepare' in sys.argv:
        dump(OUT/'plan.json',plan);dump(OUT/'warm_start_u3.json',warm_source(json.loads(SOURCE.read_text())))
        print(json.dumps(plan));return
    assert json.loads((OUT/'plan.json').read_text())==json.loads(json.dumps(plan))
    start=time.monotonic();deadline=start+240;results=[]
    for i,(seed,std) in enumerate(zip(plan['seeds'],plan['perturbation_stds'])):
        schedule=SCHEDULES[i]
        cxm=[torch.eye(16,dtype=torch.complex128) if s is None else matrix_for_schedule([s])[0] for s in schedule]
        local_start=time.monotonic();local_deadline=min(deadline,local_start+120)
        initial=x+np.random.default_rng(seed).normal(0,std,x.shape)
        best={'loss':float('inf'),'x':initial.copy(),'nfev':0}
        def fun(z):
            if time.monotonic()>=local_deadline:raise BudgetExceeded()
            loss,gradient=objective(z,fixed,cxm,v,True);best['nfev']+=1
            if loss<best['loss']:best.update(loss=loss,x=z.copy())
            return loss,gradient
        try:
            fit=minimize(fun,initial,jac=True,method='L-BFGS-B',options={'maxiter':300,'ftol':1e-15,'gtol':1e-11,'maxls':30})
            optimizer={'success':bool(fit.success),'message':str(fit.message),'nit':int(fit.nit),'nfev':int(fit.nfev),'timed_out':False}
        except BudgetExceeded:
            optimizer={'success':False,'message':'bounded compute deadline reached','nit':None,'nfev':best['nfev'],'timed_out':True}
        candidate=serialize(best['x'],schedule);path=OUT/f'case{i}.json';dump(path,candidate)
        checked=evaluate(candidate);ind=independent(candidate)
        assert abs(ind['off_diagonal_error']-checked['off_diagonal_error'])<1e-10
        results.append({'seed':seed,'std':std,'best_evaluated_loss':best['loss'],'candidate':str(path.relative_to(ROOT)),'optimizer':optimizer,'checker':checked,'independent':ind,'elapsed_seconds':time.monotonic()-local_start})
        dump(OUT/'result.json',{'plan':plan,'results':results,'elapsed_seconds':time.monotonic()-start,'complete':len(results)==2,'limitations':'two bounded starts only; failed candidates do not exclude topology or establish lower bound'})
    dump(OUT/'independent_validation.json',{'method':'separate NumPy Kronecker gates and bit permutation reconstruction from saved gate lists','results':[r['independent'] for r in results]})
    print(json.dumps({'complete':True,'elapsed_seconds':time.monotonic()-start,'results':[{k:r[k] for k in ['seed','best_evaluated_loss','independent','optimizer']} for r in results]}))

if __name__=='__main__':main()
