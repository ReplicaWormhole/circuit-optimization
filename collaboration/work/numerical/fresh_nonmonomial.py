"""One bounded full-local 13-CNOT experiment; failed fits imply no exclusion."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,hashlib,sqlite3,time,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
OUT=Path(__file__).resolve().parent
SEED=9261303
SCHEDULES=[[(0,1),(2,3),(1,2),(3,0),(0,2),(2,0),(0,2),(1,3),(3,1),(0,1),(1,0),(2,3),(3,2)],[(0,3),(1,2),(2,0),(3,1),(2,3),(0,2),(2,0),(1,3),(3,1),(1,0),(0,1),(3,2),(2,3)]]
def schedules_in(value):
    if isinstance(value,dict):
        gates=value.get('gates',[])
        if isinstance(gates,list):
            edges=[(g['control'],g['target']) for g in gates if isinstance(g,dict) and g.get('gate')=='cx']
            if len(edges)==13:yield tuple(edges)
        for key,v in value.items():
            if key in ('schedule','topology') and isinstance(v,list) and len(v)==13 and all(isinstance(x,list) and len(x)==2 and all(isinstance(y,int) for y in x) for x in v):yield tuple(map(tuple,v))
            yield from schedules_in(v)
    elif isinstance(value,list):
        for v in value:yield from schedules_in(v)
def prepare():
    prior=set(); scanned=0
    for p in ROOT.rglob('*.json'):
        if '.git' in p.parts or OUT in p.parents:continue
        try: value=json.loads(p.read_text())
        except (ValueError,OSError,UnicodeError):continue
        scanned+=1;prior.update(schedules_in(value))
    db=sqlite3.connect(ROOT/'experiments.sqlite3')
    # Config text is also searched independently of saved candidate lists.
    columns=[r[1] for r in db.execute('pragma table_info(attempts)')]
    configs=0
    if 'config_json' in columns:
        for (raw,) in db.execute('select config_json from attempts'):
            configs+=1;prior.update(schedules_in(json.loads(raw)))
    db.close()
    for s in SCHEDULES:assert tuple(s) not in prior,'Prior exact ordered schedule found'
    plan={'seed':SEED,'schedules':SCHEDULES,'rule':'Five different round-robin/nonmonomial prefix CNOTs, then eight pair-block tail CNOTs. All 14 local layers have arbitrary SU(2) gates on all four wires; no fixed target rows or prefix matrices. Two distinct fresh random starts.','maxiter_per_start':80,'wall_cap_seconds':120,'threads':1,'prior_json_files_scanned':scanned,'prior_exact_13_schedules':len(prior),'ledger_configs_scanned':configs,'novelty_scope':'Exact ordered schedules in readable saved JSON gate lists/schedule fields and ledger configuration; not equivalence classes or unsaved trials','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'initial_angles':'Independent normal(0,0.7), seed+case'}
    path=OUT/'plan.json';assert not path.exists();path.write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan))
def independent(candidate):
    import numpy as np
    u=np.eye(16,dtype=complex)
    for g in candidate['gates']:
        if g['gate']=='cx':
            m=np.zeros((16,16),complex)
            for j in range(16):
                bits=[(j>>(3-q))&1 for q in range(4)];bits[g['target']]^=bits[g['control']]
                k=sum(b<<(3-q) for q,b in enumerate(bits));m[k,j]=1
        else:
            t,p,l=[g[k] for k in ('theta','phi','lam')];c,s=np.cos(t/2),np.sin(t/2)
            a=np.array([[c,-np.exp(1j*l)*s],[np.exp(1j*p)*s,np.exp(1j*(p+l))*c]])
            m=np.array([[1]],complex)
            for q in range(4):m=np.kron(m,a if q==g['qubit'] else np.eye(2))
        u=m@u
    v=np.zeros((16,16),complex)
    for j in range(16):
        b=[(j>>(3-q))&1 for q in range(4)];b=b[-1:]+b[:-1];v[sum(x<<(3-q) for q,x in enumerate(b)),j]=1
    d=u@v@u.conj().T;off=d-np.diag(np.diag(d))
    return {'loss':float(np.sum(abs(off)**2)/16),'off_diagonal_error':float(np.max(abs(off))),'unitarity_error':float(np.max(abs(u@u.conj().T-np.eye(16))))}
def run():
    import numpy as np,torch
    from scipy.optimize import minimize
    from search13_adaptive_topology import objective,matrix_for_schedule,serialize
    from check_circuit import evaluate,right_shift
    plan=json.loads((OUT/'plan.json').read_text());assert plan['source_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    torch.set_num_threads(1);fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    rows=[];start=time.monotonic()
    class Cap(Exception):pass
    def expired(signum,frame):raise Cap()
    signal.signal(signal.SIGALRM,expired)
    for case,schedule in enumerate(plan['schedules']):
        x0=np.random.default_rng(SEED+case).normal(0,.7,(14,4,3)).ravel();best=[float('inf'),x0.copy()];calls=[0]
        cx=matrix_for_schedule(schedule)
        def f(x):
            loss,grad=objective(x,fixed,cx,v,True);calls[0]+=1
            if loss<best[0]:best[:]=[loss,x.copy()]
            return loss,grad
        remaining=120-(time.monotonic()-start)
        signal.setitimer(signal.ITIMER_REAL,max(.01,min(55,remaining)))
        try:
            result=minimize(f,x0,jac=True,method='L-BFGS-B',options={'maxiter':80,'ftol':1e-15,'gtol':1e-11,'maxls':20})
            record={'nit':int(result.nit),'message':str(result.message),'optimizer_success':bool(result.success),'wall_cap_hit':False}
        except Cap:record={'nit':None,'message':'wall cap interrupted fitting; best evaluated point retained','optimizer_success':False,'wall_cap_hit':True}
        finally:signal.setitimer(signal.ITIMER_REAL,0)
        candidate=serialize([[] for _ in range(14)],schedule,best[1].reshape(14,4,3));path=OUT/f'case{case}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
        saved=json.loads(path.read_text());a=evaluate(saved);b=independent(saved)
        assert a['cnot_count']==13 and abs(a['off_diagonal_error']-b['off_diagonal_error'])<1e-11 and abs(best[0]-b['loss'])<1e-10
        row={'case':case,'seed':SEED+case,'candidate':str(path.relative_to(ROOT)),'evaluations':calls[0],**record,'loss':b['loss'],'independent':b,'checker':a};rows.append(row);print(json.dumps(row),flush=True)
    result={'plan':'plan.json','elapsed_seconds':time.monotonic()-start,'rows':rows,'best_candidate':min(rows,key=lambda r:r['loss'])['candidate'],'scope':'Two random starts on two ordered 13-CNOT schedules with full local freedom. Numerical unsuccessful fits are not topology exclusions or a lower bound.'};(OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':prepare() if sys.argv[1]=='prepare' else run()
