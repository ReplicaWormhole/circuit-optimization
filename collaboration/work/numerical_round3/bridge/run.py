"""Bounded bridge: free13 -> commutant13 -> equivalent12 -> polish12."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy.optimize import minimize
from delete14_search import axis_rotation,u3_to_axis
from search13_adaptive_topology import objective,matrix_for_schedule
from search13_fresh_chain_dynamic_rows import decode,unitary
from collaboration.work.numerical_round3.run import warm_source,serialize,independent,dump,BudgetExceeded,matrix_u3
from check_circuit import evaluate,right_shift
OUT=Path(__file__).resolve().parent
SOURCE=ROOT/'collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json'
ORIGINAL=[(0,1),(3,2),(1,2),(1,0),(2,3),(2,1),(0,1),(3,2),(1,2),(1,0),(2,3),(2,1),(0,1)]
torch.set_num_threads(1)
def plan():
    schedules=[]
    for k in [0,6]:
        schedule=ORIGINAL.copy();schedule[k:k+3]=[(0,1),(1,2),(0,1)];schedules.append(schedule)
    return {'motif_locations':[0,6],'schedules13':schedules,'seeds':[926401,926402],'warm_noise_std':.08,'stages':['free13 max100iterations','hardcommuting13 max100iterations','collapse equivalencecheck','free12 polishmax200iterations'],'constraints':'At layersk+1 andk+2 wire0axis[x,y]=0;wire1axis[y,z]=0; all other coordinates free','collapse':'template[None,12,02]; middlelayerk+2identity;nextlayerk+3 becomes Lk+3@Lk+2','per_motif_wall_cap_seconds':120,'total_wall_cap_seconds':240,'threads':1,'optimizer':'SciPy L-BFGS-B analyticTorch gradient','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependency_hashes':{str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['collaboration/work/numerical_round3/run.py','delete14_search.py','search13_adaptive_topology.py','search13_fresh_chain_dynamic_rows.py','check_circuit.py']}}
def cxm(schedule):return [torch.eye(16,dtype=torch.complex128) if edge is None else matrix_for_schedule([edge])[0] for edge in schedule]
def fit(initial,schedule,mask,maxiter,deadline):
    initial=initial.copy();initial[~mask]=0
    best={'loss':float('inf'),'x':initial.copy(),'nfev':0}
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)]
    v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    ent=cxm(schedule)
    def fun(z):
        if time.monotonic()>=deadline:raise BudgetExceeded()
        x=np.zeros(168);x[mask]=z
        loss,g=objective(x,fixed,ent,v,True);best['nfev']+=1
        if loss<best['loss']:best.update(loss=loss,x=x.copy())
        return loss,g[mask]
    try:
        result=minimize(fun,initial[mask],jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-15,'gtol':1e-11,'maxls':30})
        best['optimizer']={'nit':int(result.nit),'success':bool(result.success),'message':str(result.message),'timeout':False}
    except BudgetExceeded:best['optimizer']={'nit':None,'success':False,'message':'deadline','timeout':True}
    if not np.isfinite(best['loss']):best['loss']=objective(initial,fixed,ent,v,False)
    return best

def save(name,x,schedule):
    candidate=serialize(x,schedule);dump(OUT/(name+'.json'),candidate)
    return {'candidate':str((OUT/(name+'.json')).relative_to(ROOT)),'checker':evaluate(candidate),'independent':independent(candidate)}
def main():
    frozen=plan()
    if '--prepare' in sys.argv:dump(OUT/'plan.json',frozen);print(json.dumps(frozen));return
    assert json.loads((OUT/'plan.json').read_text())==json.loads(json.dumps(frozen))
    source=warm_source(json.loads(SOURCE.read_text()));assert evaluate(source)['valid_diagonalizer']
    warm=decode(source).ravel();start=time.monotonic();results=[]
    for i,k in enumerate([0,6]):
        local=time.monotonic();deadline=min(start+240,local+120);schedule=[tuple(e) for e in frozen['schedules13'][i]]
        x=warm+np.random.default_rng(frozen['seeds'][i]).normal(0,.08,warm.shape)
        free=fit(x,schedule,np.ones(168,bool),100,deadline);freegate=save(f'case{i}_free13',free['x'],schedule)
        mask=np.ones((14,4,3),bool)
        for layer in [k+1,k+2]:mask[layer,0,[0,1]]=False;mask[layer,1,[1,2]]=False
        constrained=fit(free['x'],schedule,mask.ravel(),100,deadline);congate=save(f'case{i}_constrained13',constrained['x'],schedule)
        angles=constrained['x'].reshape(14,4,3).copy();old=angles.copy()
        for wire in range(4):
            a=axis_rotation(*torch.tensor(old[k+2,wire])).numpy();b=axis_rotation(*torch.tensor(old[k+3,wire])).numpy()
            angles[k+3,wire]=u3_to_axis(matrix_u3(b@a,wire))
        angles[k+2]=0;schedule12=schedule.copy();schedule12[k:k+3]=[None,(1,2),(0,2)]
        a=unitary(torch.tensor(old),cxm(schedule)).numpy();b=unitary(torch.tensor(angles),cxm(schedule12)).numpy()
        phase=np.vdot(a,b)/16;equivalence=float(np.max(abs(b-phase*a)))
        assert equivalence<1e-10
        pregate=save(f'case{i}_collapsed12',angles.ravel(),schedule12)
        polish=fit(angles.ravel(),schedule12,np.ones(168,bool),200,deadline);postgate=save(f'case{i}_final12',polish['x'],schedule12)
        records={'motif':k,'seed':frozen['seeds'][i],'free13':freegate,'constrained13':congate,'collapsed12':pregate,'final12':postgate,'equivalence_error_up_to_global_phase':equivalence,'phase_modulus':float(abs(phase)),'elapsed_seconds':time.monotonic()-local,'fits':{name:{key:value for key,value in fitresult.items() if key!='x'} for name,fitresult in [('free',free),('constrained',constrained),('polish',polish)]}}
        results.append(records);dump(OUT/'result.json',{'plan':frozen,'results':results,'elapsed_seconds':time.monotonic()-start,'limitations':'boundednumericalfits only; no exclusion or exactness'})
        print(json.dumps({'motif':k,'equivalence':equivalence,'final12':postgate['independent'],'elapsed':records['elapsed_seconds']}),flush=True)
if __name__=='__main__':main()
