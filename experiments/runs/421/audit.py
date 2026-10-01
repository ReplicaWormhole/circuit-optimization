#!/usr/bin/env python3
"""Reserved independent complete endpoint audit; no optimizer or derivative."""
from __future__ import annotations
import hashlib,inspect,json,math,os
from pathlib import Path
import signal,sqlite3,sys,time
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import scipy
from scipy.linalg import expm
import scipy.linalg._matfuncs_expm
import threadpoolctl
HERE=Path(__file__).resolve().parent
FIXED={'schema_version':1,'purpose':'independent two-endpoint audit of reordered420 ten-CX fits','method':'independent primitive and SciPy-expm reordered endpoint provenance audit','board_id':154,'board_owner':'verifier10reorder','ledger_agent':'verifier10reorder','parent_run':420,'target':'V4','source_run':420,'source_directory':'experiments/runs/420','seed':0,'fit_seed':10012301,'noise_sigmas':[.025,.25],'case_count':2,'coordinate_count':132,'cx_count':10,'matrix_build_count':14,'fit_noise_normal_draw_count':264,'time_limit_seconds':60,'thread_count':1,'matrix_tolerance':1e-9,'loss_tolerance':1e-10,'diagnostic_tolerance':1e-10,'numeric_validity_tolerance':1e-9,'evidence_scope':'Finite savedendpoint validation only; no optimizer,AD,recognition,exactcertificate,exclusion or incumbentchange.'}
REQUIRED=set(FIXED)|{'script_sha256','input_manifest_sha256','runtime_versions','dependency_sha256','plan_sha256'}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def versions():return {'python':'.'.join(map(str,sys.version_info[:3])),'numpy':np.__version__,'scipy':scipy.__version__,'threadpoolctl':threadpoolctl.__version__}
def dependencies():return {'python':sha(sys.executable),'numpy_init':sha(np.__file__),'numpy_multiarray':sha(np._core._multiarray_umath.__file__),'scipy_init':sha(scipy.__file__),'scipy_expm_source':sha(inspect.getsourcefile(expm)),'scipy_expm_core':sha(scipy.linalg._matfuncs_expm.__file__),'threadpoolctl':sha(threadpoolctl.__file__)}
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

def metrics(U,V):
    A=U@V@U.conj().T;diagonal=np.diag(A);off=A-np.diag(diagonal)
    unitary=float(np.max(np.abs(U@U.conj().T-np.eye(16))))
    maxoff=float(np.max(np.abs(off)));root=float(np.max(np.abs(diagonal**4-1)))
    return {'loss':float(np.sum(np.abs(off)**2)/16),'unitarity_error':unitary,'off_diagonal_error':maxoff,'eigenvalue_root_error':root,'valid_diagonalizer':max(unitary,maxoff,root)<=FIXED['numeric_validity_tolerance']}

def guard():
    if not HERE.name.isdigit() or HERE.parent.name!='runs' or HERE.parent.parent.name!='experiments':raise RuntimeError('reserved run directory required')
    config_raw=(HERE/'config.json').read_bytes();cfg=json.loads(config_raw)
    if config_raw!=json.dumps(cfg,sort_keys=True,separators=(',',':'),allow_nan=False).encode():raise RuntimeError('canonical raw config required')
    if sha(HERE/'PLAN.md')!=cfg['plan_sha256']:raise RuntimeError('plan hash mismatch')
    if set(cfg)!=REQUIRED:raise RuntimeError('config schema mismatch')
    for k,v in FIXED.items():
        if cfg[k]!=v or type(cfg[k]) is not type(v):raise RuntimeError('config value mismatch '+k)
    if sha(__file__)!=cfg['script_sha256'] or sha(HERE/'input_manifest.json')!=cfg['input_manifest_sha256']:raise RuntimeError('script/input manifest hash mismatch')
    if versions()!=cfg['runtime_versions'] or dependencies()!=cfg['dependency_sha256']:raise RuntimeError('runtime dependency mismatch')
    root,run=HERE.parents[2],int(HERE.name)
    if (root/'code_snapshots'/cfg['script_sha256']).read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('code snapshot mismatch')
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for rel,h in manifest['files'].items():
        path=root/rel
        if sha(path)!=h:raise RuntimeError('source artifact hash mismatch '+rel)
    with sqlite3.connect(f'file:{root/"experiments.sqlite3"}?mode=ro',uri=True) as db:
        row=db.execute('SELECT status,code_sha256,config_json,method,agent,parent_id,target FROM attempts WHERE id=?',(run,)).fetchone()
        source=db.execute('SELECT status,code_sha256,config_json,agent,parent_id,target FROM attempts WHERE id=?',(cfg['source_run'],)).fetchone()
    if row is None or row[0]!='running' or row[1]!=cfg['script_sha256'] or row[2].encode()!=config_raw or row[3]!=cfg['method'] or row[4]!=cfg['ledger_agent'] or row[5]!=cfg['parent_run'] or row[6]!=cfg['target']:raise RuntimeError('raw-config active ledger reservation mismatch')
    if source is None or source[0]!='failed' or source[1]!=manifest['fit_script_sha256'] or source[3:]!=('numerical10reorder',416,'V4'):raise RuntimeError('source run not closed with frozen script')
    source_config=json.loads((root/cfg['source_directory']/'config.json').read_text())
    if source[2].encode()!=(root/cfg['source_directory']/'config.json').read_bytes() or json.loads(source[2])!=source_config:raise RuntimeError('source ledger/config mismatch')
    with sqlite3.connect(f'file:{root/"collaboration/board.sqlite3"}?mode=ro',uri=True) as db:board=db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?',(cfg['board_id'],)).fetchone()
    if board!=('running',cfg['board_owner'],run):raise RuntimeError('board owner/run link mismatch')
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        if os.environ[key]!='1':raise RuntimeError('thread environment mismatch')
    return root,run,cfg,manifest

def checked_vector(raw):
    out=np.asarray(raw,dtype=float)
    if out.shape!=(132,) or not np.isfinite(out).all():raise AssertionError('vector is not finite full132')
    return out

def check_word(candidate,topology):
    if candidate['n']!=4 or len(candidate['gates'])!=54:raise AssertionError('endpoint must have54gates')
    expected=[]
    for k in range(11):
        expected.extend(('u3',q) for q in range(4))
        if k<10:expected.append(('cx',*topology[k]))
    actual=[('cx',g['control'],g['target']) if g['gate']=='cx' else (g['gate'],g['qubit']) for g in candidate['gates']]
    if actual!=expected:raise AssertionError('complete endpoint chronology/topology differs')

def main():
    started=time.monotonic();root,run,cfg,manifest=guard()
    with (HERE/'execution_started.json').open('x') as marker:json.dump({'run_id':run,'config_sha256':sha(HERE/'config.json')},marker)
    def stop(signum,frame):raise TimeoutError('60-second endpoint audit budget')
    signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,max(.001,60-(time.monotonic()-started)))
    records=[];builds=0
    try:
        limiter=threadpoolctl.threadpool_limits(limits=1);src=root/cfg['source_directory'];load=lambda rel:json.loads((src/rel).read_text())
        result=load('result.json');fitcfg=load('config.json');source=load('source_case_result.json');oldword=load('source_best_candidate.json');plans=load('topologies.json');preflight=json.loads((root/'experiments/runs/419/result.json').read_text())
        if not preflight['passed'] or preflight['matrix_build_count']!=22 or preflight['point_count']!=3:raise AssertionError('successful419 preflight required')
        if result['run_id']!=420 or result['status']!='completed' or result['case_count']!=2 or len(result['records'])!=2:raise AssertionError('source completion/order mismatch')
        for k,expected in [('seed',cfg['fit_seed']),('noise_sigmas',cfg['noise_sigmas']),('case_count',2),('free_parameter_count',132),('expected_cx_count',10),('local_layer_count',11),('maxiter',250),('maxfun',10000),('parent_run',416),('target','V4')]:
            if fitcfg[k]!=expected:raise AssertionError('fitconfig mismatch '+k)
        if result['script_sha256']!=manifest['fit_script_sha256'] or result['config_sha256']!=manifest['fit_config_sha256'] or result['family_sha256']!=manifest['family_sha256'] or result['source_case_result_sha256']!=manifest['source_case_result_sha256']:raise AssertionError('source result hashbindings differ')
        topology=plans['new_topology'];oldtopology=plans['old_topology']
        if topology!=fitcfg['new_topology'] or oldtopology!=fitcfg['old_topology'] or source['delete_cx_ordinal']!=7 or source['topology']!=oldtopology:raise AssertionError('source/new topology mismatch')
        if load('base_vector.json')['base132']!=source['best132']:raise AssertionError('rawbase altered')
        literal=json.loads(json.dumps(oldword));positions=[j for j,g in enumerate(literal['gates']) if g['gate']=='cx']
        if [[literal['gates'][j]['control'],literal['gates'][j]['target']] for j in positions]!=oldtopology:raise AssertionError('source oldword edge mismatch')
        g6,g7=literal['gates'][positions[6]],literal['gates'][positions[7]];g6['control'],g7['control']=g7['control'],g6['control'];g6['target'],g7['target']=g7['target'],g6['target']
        check_word(literal,topology);newbase=load('base_candidate.json')
        if literal['n']!=newbase['n'] or literal['gates']!=newbase['gates']:raise AssertionError('semantic base transplant mismatch')
        rng=np.random.Generator(np.random.PCG64(cfg['fit_seed']));noise=[rng.normal(0,sigma,size=132) for sigma in cfg['noise_sigmas']]
        if [n.tolist() for n in noise]!=load('all_noise132.json') or result['noise_sha256']!=sha(src/'all_noise132.json'):raise AssertionError('noiseallocation/hash mismatch')
        V=shift()
        for ordinal,sigma in enumerate(cfg['noise_sigmas']):
            path=Path('cases')/f'{ordinal:02d}_sigma{sigma:g}';case=load(path/'result.json');vectors=load(path/'vectors.json')
            if case!=result['records'][ordinal] or case['case_index']!=ordinal or case['noise_sigma']!=sigma or case['topology']!=topology or case['old_topology']!=oldtopology:raise AssertionError('caseorder/topology mismatch')
            base=checked_vector(vectors['base132']);initial=checked_vector(vectors['initial132']);best=checked_vector(case['best132']);savednoise=checked_vector(vectors['noise132'])
            if vectors['base132']!=source['best132'] or vectors['noise_sigma']!=sigma or vectors['new_topology']!=topology or vectors['old_topology']!=oldtopology or not np.array_equal(savednoise,noise[ordinal]) or not np.array_equal(initial,base+noise[ordinal]):raise AssertionError('rawbase/noise/initial mismatch')
            if vectors['base132']!=preflight['records'][0]['vector132'] or vectors['initial132']!=preflight['records'][ordinal+1]['vector132'] or vectors['noise132']!=preflight['records'][ordinal+1]['noise132']:raise AssertionError('419coordinates/noise mismatch')
            savedliteral=load(path/'literal_base_candidate.json')
            if literal['n']!=savedliteral['n'] or literal['gates']!=savedliteral['gates']:raise AssertionError('case literalbase changed')
            words=[load(path/name) for name in ['base_candidate.json','initial_candidate.json','best_candidate.json']]
            for word in words:check_word(word,topology)
            L=gate_word(literal);vector_matrices=[vector_word(v,topology) for v in [base,initial,best]];serialized_matrices=[gate_word(word) for word in words];builds+=7
            errors={'literal_to_base':phase_error(L,vector_matrices[0]),'base_serializer':phase_error(vector_matrices[0],serialized_matrices[0]),'initial_serializer':phase_error(vector_matrices[1],serialized_matrices[1]),'best_serializer':phase_error(vector_matrices[2],serialized_matrices[2])}
            diagnostics=[metrics(U,V) for U in vector_matrices];loss_errors={};metric_errors={}
            for label,ours,saved in [('base',diagnostics[0],case['base_checker']),('initial',diagnostics[1],case['initial_checker']),('best',diagnostics[2],case['best_checker'])]:
                loss_errors[label]=abs(ours['loss']-case[label+'_loss'])
                for key in ['unitarity_error','off_diagonal_error','eigenvalue_root_error']:metric_errors[label+'_'+key]=abs(ours[key]-saved[key])
                if ours['valid_diagonalizer']!=saved['valid_diagonalizer'] or saved['cnot_count']!=10:raise AssertionError('validity/count differs')
            if load(path/'base_metrics.json')!={'loss':case['base_loss'],'checker':case['base_checker']} or load(path/'initial_metrics.json')!={'loss':case['initial_loss'],'checker':case['initial_checker']}:raise AssertionError('metricsnapshot differs')
            if case['candidate_file']!=str(path/'best_candidate.json') or case['candidate_sha256']!=sha(src/case['candidate_file']):raise AssertionError('candidatefile/hash differs')
            history=case['history']
            if len(history)!=case['evaluations'] or case['evaluations']>10000 or case['iterations']>250 or case['free_parameter_count']!=132 or [h['evaluation'] for h in history]!=list(range(1,len(history)+1)):raise AssertionError('history/countbudget differs')
            if abs(case['best_loss']-min([case['initial_loss']]+[h['loss'] for h in history]))>cfg['loss_tolerance']:raise AssertionError('bestlosshistory differs')
            passed=max(errors.values())<=cfg['matrix_tolerance'] and max(loss_errors.values())<=cfg['loss_tolerance'] and max(metric_errors.values())<=cfg['diagnostic_tolerance'] and max(x['unitarity_error'] for x in diagnostics)<=cfg['matrix_tolerance']
            records.append({'case_index':ordinal,'noise_sigma':sigma,'matrix_phase_errors':errors,'loss_errors':loss_errors,'diagnostic_errors':metric_errors,'base_initial_best_metrics':diagnostics,'iterations':case['iterations'],'evaluations':case['evaluations'],'termination':case['termination'],'passed':bool(passed)});dump(HERE/'partial_result.json',{'matrix_build_count':builds,'records':records})
        best=min(result['records'],key=lambda r:(r['best_loss'],r['case_index']));minmaxoff=min(result['records'],key=lambda r:(r['best_checker']['off_diagonal_error'],r['case_index']))
        if result['best']!=best or load('best_candidate.json')!=load(best['candidate_file']) or result['best_candidate_sha256']!=sha(src/'best_candidate.json') or result['adaptive_ranking'] is not False or result['ranking']!=fitcfg['ranking']:raise AssertionError('deterministicglobalminimumsummary differs')
        if builds!=14:raise AssertionError('finite matrixcount differs')
        output={'run_id':run,'source_run':420,'passed':all(r['passed'] for r in records),'case_count':len(records),'matrix_build_count':builds,'fit_noise_normal_draw_count':264,'minimum_loss_case_index':best['case_index'],'minimum_maxoff_case_index':minmaxoff['case_index'],'minimum_loss_differs_minimum_maxoff':best['case_index']!=minmaxoff['case_index'],'both_best_endpoints_invalid':all(not r['base_initial_best_metrics'][2]['valid_diagonalizer'] for r in records),'source_total_iterations':sum(r['iterations'] for r in result['records']),'source_total_evaluations':sum(r['evaluations'] for r in result['records']),'script_sha256':cfg['script_sha256'],'config_sha256':sha(HERE/'config.json'),'input_manifest_sha256':cfg['input_manifest_sha256'],'runtime_versions':versions(),'dependency_sha256':dependencies(),'records':records,'seconds':time.monotonic()-started,'evidence_scope':cfg['evidence_scope']};dump(HERE/'result.json',output);print(json.dumps({k:output[k] for k in ['run_id','passed','matrix_build_count','case_count','minimum_loss_case_index','minimum_maxoff_case_index','seconds']}));return 0 if output['passed'] else 1
    except BaseException as error:
        dump(HERE/'failure.json',{'run_id':run,'error_type':type(error).__name__,'error':str(error),'completed_cases':len(records),'matrix_build_count':builds,'seconds':time.monotonic()-started});raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        if 'limiter' in locals():limiter.restore_original_limits()
if __name__=='__main__':raise SystemExit(main())
