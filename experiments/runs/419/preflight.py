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
FIXED={'schema_version':1,'purpose':'independent finite preflight of transplanted base and two reordered ten-CX starts','method':'independent primitive and SciPy-expm reordered-word directional autograd audit','board_id':153,'board_owner':'verifier10reorder','ledger_agent':'verifier10reorder','parent_run':416,'target':'V4','seed':0,'fit_seed':10012301,'noise_sigmas':[.025,.25],'start_count':2,'point_count':3,'coordinate_count':132,'cx_count':10,'matrix_build_count':22,'autograd_gradient_count':3,'finite_difference_loss_count':6,'fit_noise_normal_draw_count':264,'direction_random_draw_count':0,'time_limit_seconds':60,'thread_count':1,'matrix_tolerance':1e-9,'loss_tolerance':1e-10,'gradient_absolute_tolerance':1e-6,'gradient_relative_tolerance':1e-5,'finite_difference_step':1e-6,'topology':[[0,2],[1,3],[0,1],[2,3],[0,2],[3,2],[1,0],[0,2],[1,2],[3,2]],'direction_rule':'d_j=cos((j+1)*(point_ordinal+1)/37),j0..131; Euclidean normalize; no RNG','evidence_scope':'Finite three-point implementation/initialization preflight only; no fit, exact proof, exclusion or incumbent change.'}
REQUIRED=set(FIXED)|{'script_sha256','input_manifest_sha256','runtime_versions','dependency_sha256','plan_sha256'}
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
    raw=(HERE/'config.json').read_bytes();cfg=json.loads(raw)
    if raw!=json.dumps(cfg,sort_keys=True,separators=(',',':'),allow_nan=False).encode() or set(cfg)!=REQUIRED:raise RuntimeError('canonical raw config schema required')
    for k,v in FIXED.items():
        if cfg[k]!=v or type(cfg[k]) is not type(v):raise RuntimeError('fixed config mismatch '+k)
    if sha(__file__)!=cfg['script_sha256'] or sha(HERE/'PLAN.md')!=cfg['plan_sha256'] or sha(HERE/'input_manifest.json')!=cfg['input_manifest_sha256']:raise RuntimeError('source/plan/manifest mismatch')
    if versions()!=cfg['runtime_versions'] or dependencies()!=cfg['dependency_sha256']:raise RuntimeError('runtime dependency mismatch')
    root,run=HERE.parents[2],int(HERE.name);manifest=json.loads((HERE/'input_manifest.json').read_text())
    if (root/'code_snapshots'/cfg['script_sha256']).read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('reserved code snapshot mismatch')
    for rel,h in manifest['files'].items():
        if sha(root/rel)!=h:raise RuntimeError('input hash mismatch '+rel)
    for name,key in [('fit_family.py','family_sha256'),('fit_source.py','fit_script_sha256'),('fit_config.json','fit_config_sha256')]:
        if sha(HERE/name)!=manifest[key]:raise RuntimeError('copied audit-subject hash mismatch '+name)
    with sqlite3.connect(f'file:{root/"experiments.sqlite3"}?mode=ro',uri=True) as db:
        row=db.execute('SELECT status,code_sha256,config_json,method,agent,parent_id,target FROM attempts WHERE id=?',(run,)).fetchone()
        source=db.execute('SELECT status,agent,parent_id,target FROM attempts WHERE id=416').fetchone()
    if row!=('running',cfg['script_sha256'],raw.decode(),cfg['method'],cfg['ledger_agent'],cfg['parent_run'],cfg['target']):raise RuntimeError('raw-config active ledger reservation mismatch')
    if source!=('failed','numerical10',411,'V4'):raise RuntimeError('closed source416 provenance required')
    with sqlite3.connect(f'file:{root/"collaboration/board.sqlite3"}?mode=ro',uri=True) as db:board=db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?',(cfg['board_id'],)).fetchone()
    if board!=('running',cfg['board_owner'],run):raise RuntimeError('matching board ownership/link required')
    return root,run,cfg,manifest

def main():
    started=time.monotonic();root,run,cfg,manifest=guard()
    with (HERE/'execution_started.json').open('x') as marker:json.dump({'run_id':run,'config_sha256':sha(HERE/'config.json')},marker)
    def stop(signum,frame):raise TimeoutError('60-second preflight budget')
    signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,max(.001,60-(time.monotonic()-started)))
    records=[];builds=0;gradients=0
    try:
        limiter=threadpoolctl.threadpool_limits(limits=1);torch.set_num_threads(1);torch.set_num_interop_threads(1)
        source=json.loads((root/'experiments/runs/416/cases/07_delete_cx07/result.json').read_text())
        fitcfg=json.loads((HERE/'fit_config.json').read_text())
        for key,expected in [('seed',cfg['fit_seed']),('noise_sigmas',cfg['noise_sigmas']),('case_count',2),('free_parameter_count',132),('local_layer_count',11),('expected_cx_count',10),('new_topology',cfg['topology']),('parent_run',416),('target','V4')]:
            if fitcfg[key]!=expected:raise AssertionError('fit config mismatch '+key)
        savedbase=json.loads((root/'collaboration/work/numerical10reorder/base_vector.json').read_text())
        if savedbase['base132']!=source['best132'] or source['delete_cx_ordinal']!=7:raise AssertionError('raw source base altered')
        base=np.asarray(source['best132'],float);topology=cfg['topology']
        if base.shape!=(132,) or not np.isfinite(base).all():raise AssertionError('raw full132 base required')
        saved=json.loads((root/'experiments/runs/416/cases/07_delete_cx07/best_candidate.json').read_text())
        literal=json.loads(json.dumps(saved));positions=[j for j,g in enumerate(literal['gates']) if g['gate']=='cx'];oldedges=[[saved['gates'][j]['control'],saved['gates'][j]['target']] for j in positions]
        if oldedges!=source['topology'] or oldedges!=fitcfg['old_topology']:raise AssertionError('source old chronology mismatch')
        g6,g7=literal['gates'][positions[6]],literal['gates'][positions[7]]
        g6['control'],g7['control']=g7['control'],g6['control'];g6['target'],g7['target']=g7['target'],g6['target']
        if [[literal['gates'][j]['control'],literal['gates'][j]['target']] for j in positions]!=topology:raise AssertionError('transplant edge swap differs')
        newbase=json.loads((root/'collaboration/work/numerical10reorder/base_candidate.json').read_text())
        if literal['n']!=newbase['n'] or literal['gates']!=newbase['gates']:raise AssertionError('literal base transplant gate provenance mismatch')
        if len(literal['gates'])!=54 or sum(g['gate']=='cx' for g in literal['gates'])!=10:raise AssertionError('transplant word count mismatch')
        L=gate_word(literal);builds+=1
        rng=np.random.Generator(np.random.PCG64(cfg['fit_seed']));noise=[rng.normal(0,sigma,size=132) for sigma in cfg['noise_sigmas']]
        points=[base,base+noise[0],base+noise[1]];family=load_family();word=family.TorchWord(topology);V=shift()
        for ordinal,x in enumerate(points):
            R=vector_word(x,topology);serialized=family.serialize(x,topology);S=gate_word(serialized);N=family.unitary_numpy(x,topology);T=word.unitary(torch.tensor(x,dtype=torch.float64)).detach().numpy();builds+=4
            expected=[]
            for k in range(11):
                expected.extend(('u3',q) for q in range(4))
                if k<10:expected.append(('cx',*topology[k]))
            actual=[('cx',g['control'],g['target']) if g['gate']=='cx' else (g['gate'],g['qubit']) for g in serialized['gates']]
            if serialized['n']!=4 or actual!=expected:raise AssertionError('fullword chronology mismatch')
            errors={'serializer':phase_error(R,S),'family_numpy':phase_error(R,N),'torch':phase_error(R,T)}
            if ordinal==0:errors['transplanted_literal']=phase_error(R,L)
            direction=np.cos((np.arange(132)+1)*(ordinal+1)/37);direction/=np.linalg.norm(direction)
            tx=torch.tensor(x,dtype=torch.float64,requires_grad=True);value=word.loss(tx);gradient=torch.autograd.grad(value,tx)[0].detach().numpy();builds+=1;gradients+=1
            h=cfg['finite_difference_step'];plus=loss(vector_word(x+h*direction,topology),V);minus=loss(vector_word(x-h*direction,topology),V);builds+=2
            fd=(plus-minus)/(2*h);ad=float(np.dot(gradient,direction));gerror=abs(ad-fd);threshold=cfg['gradient_absolute_tolerance']+cfg['gradient_relative_tolerance']*max(abs(ad),abs(fd));lerror=abs(float(value.detach())-loss(R,V))
            passed=max(errors.values())<=cfg['matrix_tolerance'] and lerror<=cfg['loss_tolerance'] and gerror<=threshold and np.isfinite(gradient).all()
            record={'point_ordinal':ordinal,'point_kind':'base' if ordinal==0 else 'initial','vector132':x.tolist(),'noise132':None if ordinal==0 else noise[ordinal-1].tolist(),'complete_candidate':serialized,'direction':direction.tolist(),'matrix_phase_errors':errors,'loss':loss(R,V),'loss_error':lerror,'gradient_ad_directional':ad,'gradient_fd_directional':fd,'gradient_error':gerror,'gradient_threshold':threshold,'passed':bool(passed)};records.append(record);dump(HERE/'partial_result.json',{'records':records,'matrix_build_count':builds,'autograd_gradient_count':gradients})
        if builds!=22 or gradients!=3:raise AssertionError('finite count mismatch')
        result={'run_id':run,'passed':all(r['passed'] for r in records),'point_count':len(records),'matrix_build_count':builds,'autograd_gradient_count':gradients,'finite_difference_loss_count':6,'fit_noise_normal_draw_count':264,'script_sha256':cfg['script_sha256'],'config_sha256':sha(HERE/'config.json'),'input_manifest_sha256':cfg['input_manifest_sha256'],'runtime_versions':versions(),'dependency_sha256':dependencies(),'records':records,'transplanted_literal':literal,'seconds':time.monotonic()-started,'evidence_scope':cfg['evidence_scope']};dump(HERE/'result.json',result);print(json.dumps({k:result[k] for k in ['run_id','passed','point_count','matrix_build_count','autograd_gradient_count','seconds']}));return 0 if result['passed'] else 1
    except BaseException as error:
        dump(HERE/'failure.json',{'run_id':run,'error_type':type(error).__name__,'error':str(error),'completed_points':len(records),'matrix_build_count':builds,'seconds':time.monotonic()-started});raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        if 'limiter' in locals():limiter.restore_original_limits()
if __name__=='__main__':raise SystemExit(main())
