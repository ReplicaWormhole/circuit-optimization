#!/usr/bin/env python3
"""Reserved finite independent preflight, with no optimizer or fit."""
from __future__ import annotations
import hashlib, importlib.util, inspect, json, math, os
from pathlib import Path
import signal, sqlite3, sys, time
from types import SimpleNamespace
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import scipy
from scipy.linalg import expm
import scipy.linalg._matfuncs_expm
import torch
import threadpoolctl

HERE=Path(__file__).resolve().parent
FIXED={'schema_version':1,'purpose':'independent finite preflight of eleven full-local ten-CX deletion words','method':'independent SciPy-expm and primitive matrices plus directional autograd audit','board_id':145,'seed':0,'fit_seed':10012201,'noise_sigma':0.025,'case_count':11,'coordinate_count':132,'cx_count':10,'matrix_build_count':121,'fit_noise_normal_draw_count':1452,'direction_random_draw_count':0,'autograd_gradient_count':11,'finite_difference_loss_count':22,'time_limit_seconds':60,'thread_count':1,'matrix_tolerance':1e-9,'loss_tolerance':1e-10,'gradient_absolute_tolerance':1e-6,'gradient_relative_tolerance':1e-5,'finite_difference_step':1e-6,'direction_rule':'d_j=cos((j+1)*(deletion_ordinal+1)/37), j=0..131; normalize Euclidean norm; no random direction draws','evidence_scope':'Finite preflight consistency only; no optimization, diagonalizer, exclusion or exact proof claim.'}
FIXED.update(board_id=145, board_owner='verifier10', ledger_agent='verifier10', parent_run=411, target='V4')
REQUIRED=set(FIXED)|{'script_sha256','family_sha256','fit_script_sha256','fit_config_sha256','candidate_sha256','runtime_versions','dependency_sha256','plan_sha256','acceptance_sha256'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def versions():return {'python':'.'.join(map(str,sys.version_info[:3])),'numpy':np.__version__,'scipy':scipy.__version__,'torch':torch.__version__,'threadpoolctl':threadpoolctl.__version__}
def dependencies():return {'python':sha(sys.executable),'numpy_init':sha(np.__file__),'numpy_multiarray':sha(np._core._multiarray_umath.__file__),'scipy_init':sha(scipy.__file__),'scipy_expm_source':sha(inspect.getsourcefile(expm)),'scipy_expm_core':sha(scipy.linalg._matfuncs_expm.__file__),'torch_init':sha(torch.__file__),'torch_C':sha(torch._C.__file__),'torch_cpu':sha(Path(torch.__file__).parent/'lib/libtorch_cpu.so'),'threadpoolctl':sha(threadpoolctl.__file__)}
def angle(v):
    if isinstance(v,(int,float)):return float(v)
    raw=v.replace(' ','').replace('*','')
    numerator,_,denominator=raw.partition('/')
    factor=numerator.replace('pi','')
    return float(-1 if factor=='-' else 1 if factor in ('','+') else factor)*math.pi/int(denominator or 1)

def local(g):
    name=g['gate']
    if name=='h':return np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
    if name=='x':return np.array([[0,1],[1,0]],complex)
    if name in ('s','sdg'):return np.diag([1,1j if name=='s' else -1j])
    t=angle(g['theta']);c,s=math.cos(t/2),math.sin(t/2)
    if name=='rx':return np.array([[c,-1j*s],[-1j*s,c]],complex)
    if name=='ry':return np.array([[c,-s],[s,c]],complex)
    if name=='rz':return np.diag([np.exp(-.5j*t),np.exp(.5j*t)])
    if name=='u3':
        p,l=angle(g['phi']),angle(g['lam'])
        return np.array([[c,-np.exp(1j*l)*s],[np.exp(1j*p)*s,np.exp(1j*(p+l))*c]],complex)
    raise ValueError(name)

def cx(c,t):
    M=np.zeros((16,16),complex)
    for col in range(16):
        bits=[(col>>(3-q))&1 for q in range(4)]
        bits[t]^=bits[c]
        row=sum(b<<(3-q) for q,b in enumerate(bits));M[row,col]=1
    return M

def lift(M,q):
    out=np.ones((1,1),complex)
    for wire in range(4):out=np.kron(out,M if wire==q else np.eye(2))
    return out

def gate_word(candidate):
    U=np.eye(16,dtype=complex)
    for g in candidate['gates']:U=(cx(g['control'],g['target']) if g['gate']=='cx' else lift(local(g),g['qubit']))@U
    return U

def vector_word(v,topology):
    # Independent exponential implementation; no frozen-family rotvec helper.
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1])
    U=np.eye(16,dtype=complex)
    for k,layer in enumerate(np.asarray(v).reshape(11,4,3)):
        M=np.ones((1,1),complex)
        for r in layer:M=np.kron(M,expm(-.5j*(r[0]*X+r[1]*Y+r[2]*Z)))
        U=M@U
        if k<10:U=cx(*topology[k])@U
    return U

def shift():
    M=np.zeros((16,16),complex)
    for col in range(16):
        bits=[(col>>(3-q))&1 for q in range(4)];new=bits[-1:]+bits[:-1]
        M[sum(b<<(3-q) for q,b in enumerate(new)),col]=1
    return M

def loss(U,V):
    A=U@V@U.conj().T;np.fill_diagonal(A,0)
    return float(np.sum(A.real*A.real+A.imag*A.imag)/16)

def phase_error(U,V):
    overlap=np.trace(U@V.conj().T);phase=overlap/abs(overlap) if abs(overlap)>1e-14 else 1
    return float(np.max(np.abs(U-phase*V)))

def load_family():
    spec=importlib.util.spec_from_file_location('frozen_fit_family',HERE/'fit_family.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def guard():
    if not HERE.name.isdigit() or HERE.parent.name!='runs' or HERE.parent.parent.name!='experiments':raise RuntimeError('reserved run workspace required')
    config_raw=(HERE/'config.json').read_bytes();cfg=json.loads(config_raw)
    if config_raw!=json.dumps(cfg,sort_keys=True,separators=(',',':'),allow_nan=False).encode():raise RuntimeError('canonical raw config bytes required')
    if set(cfg)!=REQUIRED:raise RuntimeError('config key mismatch')
    for k,v in FIXED.items():
        if cfg[k]!=v or type(cfg[k]) is not type(v):raise RuntimeError('frozen config mismatch '+k)
    for name,key in [('preflight.py','script_sha256'),('fit_family.py','family_sha256'),('fit_source.py','fit_script_sha256'),('fit_config.json','fit_config_sha256'),('source_candidate.json','candidate_sha256')]:
        if sha(HERE/name)!=cfg[key]:raise RuntimeError('frozen input hash mismatch '+name)
    if sha(HERE/'PLAN.md')!=cfg['plan_sha256']:raise RuntimeError('plan hash mismatch')
    if cfg['runtime_versions']!=versions() or cfg['dependency_sha256']!=dependencies():raise RuntimeError('runtime dependency mismatch')
    root,run=HERE.parents[2],int(HERE.name)
    if cfg['acceptance_sha256']!=sha(root/'collaboration/EXACT11_ACCEPTANCE.json') or cfg['candidate_sha256']!=sha(root/'experiments/runs/411/candidate.json'):raise RuntimeError('accepted source changed')
    if (root/'code_snapshots'/cfg['script_sha256']).read_bytes()!=(HERE/'preflight.py').read_bytes():raise RuntimeError('reserved code snapshot mismatch')
    fit_cfg=json.loads((HERE/'fit_config.json').read_text())
    for field,expected in [('family_sha256',cfg['family_sha256']),('script_sha256',cfg['fit_script_sha256']),('source_candidate_sha256',cfg['candidate_sha256']),('seed',cfg['fit_seed']),('noise_sigma',cfg['noise_sigma']),('case_count',cfg['case_count']),('free_parameter_count',cfg['coordinate_count']),('expected_cx_count',cfg['cx_count']),('local_layer_count',11)]:
        if fit_cfg[field]!=expected:raise RuntimeError('fit-config consistency mismatch '+field)
    with sqlite3.connect(f'file:{root/"experiments.sqlite3"}?mode=ro',uri=True) as db:row=db.execute('SELECT status,code_sha256,config_json,method,agent,parent_id,target FROM attempts WHERE id=?',(run,)).fetchone()
    if row is None or row[0]!='running' or row[1]!=cfg['script_sha256'] or row[2].encode()!=config_raw or row[3]!=cfg['method'] or row[4]!=cfg['ledger_agent'] or row[5]!=cfg['parent_run'] or row[6]!=cfg['target']:raise RuntimeError('matching raw-config active ledger reservation required')
    with sqlite3.connect(f'file:{root/"collaboration/board.sqlite3"}?mode=ro',uri=True) as db:row=db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?',(cfg['board_id'],)).fetchone()
    if row!=('running',cfg['board_owner'],run):raise RuntimeError('matching active board link required')
    if any(os.environ[key]!='1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')):raise RuntimeError('single-thread environment mismatch')
    return run,cfg

def main():
    started=time.monotonic();run,cfg=guard()
    with (HERE/'execution_started.json').open('x') as marker:json.dump({'run_id':run,'config_sha256':sha(HERE/'config.json')},marker)
    def stop(signum,frame):raise TimeoutError('60-second preflight budget')
    signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,max(.001,60-(time.monotonic()-started)))
    records=[];builds=0;gradients=0
    try:
        limiter=threadpoolctl.threadpool_limits(limits=cfg['thread_count'])
        torch.set_num_threads(1);torch.set_num_interop_threads(1)
        family=load_family();source=json.loads((HERE/'source_candidate.json').read_text());V=shift()
        rng=np.random.Generator(np.random.PCG64(cfg['fit_seed']));noise=np.array([rng.normal(0,cfg['noise_sigma'],size=132) for ordinal in range(11)])
        for ordinal in range(11):
            base,topology,deleted,layers,position=family.collapse(source,ordinal,SimpleNamespace(one_qubit_matrix=local))
            cx_positions=[p for p,g in enumerate(source['gates']) if g['gate']=='cx']
            if position!=cx_positions[ordinal]:raise AssertionError('wrong chronological deletion ordinal')
            expected=[dict(g) for p,g in enumerate(source['gates']) if p!=position]
            if topology!=[[g['control'],g['target']] for g in expected if g['gate']=='cx']:raise AssertionError('incorrect collapsed topology')
            if deleted!={'n':4,'gates':expected} or len(topology)!=10 or np.shape(base)!=(132,):raise AssertionError('collapse literal topology mismatch')
            initial=base+noise[ordinal]
            D=gate_word(deleted);B=vector_word(base,topology);BS=gate_word(family.serialize(base,topology));BN=family.unitary_numpy(base,topology)
            I0=vector_word(initial,topology);IS=gate_word(family.serialize(initial,topology));IN=family.unitary_numpy(initial,topology)
            word=family.TorchWord(topology);IT=word.unitary(torch.tensor(initial,dtype=torch.float64)).detach().numpy();builds+=8
            errors={'deleted_to_base':phase_error(D,B),'base_serializer':phase_error(B,BS),'base_family_numpy':phase_error(B,BN),'initial_serializer':phase_error(I0,IS),'initial_family_numpy':phase_error(I0,IN),'initial_torch':phase_error(I0,IT)}
            direction=np.cos((np.arange(132)+1)*(ordinal+1)/37);direction/=np.linalg.norm(direction)
            x=torch.tensor(initial,dtype=torch.float64,requires_grad=True);value=word.loss(x);gradient=torch.autograd.grad(value,x)[0].detach().numpy();builds+=1;gradients+=1
            h=cfg['finite_difference_step'];plus=loss(vector_word(initial+h*direction,topology),V);minus=loss(vector_word(initial-h*direction,topology),V);builds+=2
            fd=(plus-minus)/(2*h);ad=float(np.dot(gradient,direction));gerror=abs(ad-fd);gthreshold=cfg['gradient_absolute_tolerance']+cfg['gradient_relative_tolerance']*max(abs(ad),abs(fd));lerror=abs(float(value.detach())-loss(I0,V))
            counts=[]
            for vector in (base,initial):
                saved=family.serialize(vector,topology);counts.append((len(saved['gates']),sum(g['gate']=='cx' for g in saved['gates'])))
            passed=max(errors.values())<=cfg['matrix_tolerance'] and lerror<=cfg['loss_tolerance'] and gerror<=gthreshold and np.isfinite(gradient).all() and counts==[(54,10),(54,10)]
            record={'deletion_ordinal':ordinal,'deleted_source_position':position,'topology':topology,'base_vector':base.tolist(),'noise_vector':noise[ordinal].tolist(),'initial_vector':initial.tolist(),'direction':direction.tolist(),'matrix_phase_errors':errors,'gradient_ad_directional':ad,'gradient_fd_directional':fd,'gradient_error':gerror,'gradient_threshold':gthreshold,'loss_error':lerror,'base_initial_gate_counts':counts,'passed':bool(passed)}
            records.append(record);dump(HERE/'partial_result.json',{'run_id':run,'records':records,'matrix_builds':builds,'gradient_evaluations':gradients})
        result={'run_id':run,'passed':all(r['passed'] for r in records),'case_count':len(records),'matrix_build_count':builds,'autograd_gradient_count':gradients,'finite_difference_loss_count':22,'source_sha256':cfg['candidate_sha256'],'script_sha256':cfg['script_sha256'],'config_sha256':sha(HERE/'config.json'),'runtime_versions':versions(),'dependency_sha256':dependencies(),'records':records,'seconds':time.monotonic()-started,'evidence_scope':cfg['evidence_scope']}
        if builds!=121 or gradients!=11:raise AssertionError('frozen count mismatch')
        dump(HERE/'result.json',result);print(json.dumps({k:result[k] for k in ('run_id','passed','case_count','matrix_build_count','autograd_gradient_count','seconds')}));return 0 if result['passed'] else 1
    except BaseException as e:
        dump(HERE/'failure.json',{'run_id':run,'error_type':type(e).__name__,'error':str(e),'completed_cases':len(records),'matrix_builds':builds,'seconds':time.monotonic()-started});raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        if 'limiter' in locals():limiter.restore_original_limits()
if __name__=='__main__':raise SystemExit(main())
