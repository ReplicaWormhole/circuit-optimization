#!/usr/bin/env python3
"""Two frozen reduced-v0 diagnostics and fullV4 comparators, requiring coordinator reservation."""
from __future__ import annotations
import os
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
for key in THREAD_KEYS:
    os.environ[key] = '1'
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import signal
import sqlite3
import sys
import time
import numpy as np
import scipy
from scipy.optimize import minimize
import scipy.optimize._lbfgsb
import torch
import threadpoolctl
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
FIXED = {'schema_version': 1,
 'purpose': 'two bounded reduced v0 reordered-four-CX diagnostics with full ten-CX V4 comparators',
 'method': 'SciPy L-BFGS-B with Torch analytic gradients; reduced Wv0 full30 x-y SU2 locals; two '
           'perturbations and fullV4 comparators',
 'board_id': 156,
 'agent': 'numerical_v0',
 'parent_run': 411,
 'target': 'V4',
 'objective_target': 'Wv0=CZxy CXu-to-x with CX acting first; not V3',
 'source_filename': 'old_source.json',
 'source_sha256': '93030d184db2725c2b32e68b8b28eb8481727abbdf9f38af9018771b32842b82',
 'prefix_filename': 'prefix.json',
 'seed': 10012401,
 'rng': 'NumPy Generator PCG64',
 'noise_sigmas': [0.025, 0.25],
 'noise_rule': 'one PCG64(seed) generator; two consecutive normal(0,sigma,size=30) draws in case '
               'order0..1 with sigma0.025 then0.25; allocate all upfront',
 'case_count': 2,
 'free_parameter_count': 30,
 'local_layer_count': 5,
 'local_wire_order': [0, 2],
 'fixed_identity_wire': 1,
 'reduced_n': 3,
 'full_n': 4,
 'expected_reduced_cx_count': 4,
 'expected_full_cx_count': 10,
 'expected_reduced_gate_count': 14,
 'expected_full_gate_count': 22,
 'old_topology': [[0, 2], [0, 2], [1, 0], [1, 2]],
 'new_topology': [[0, 2], [1, 0], [0, 2], [1, 2]],
 'topology_rule': 'swap old CX ordinals1 and2; transplant all five x/y local layers unchanged; u locals '
                  'remain identity',
 'full_rule': 'accepted411 first6 Bell/router gates, then new4CX x-y local layers with bareCX32 '
              'immediately after firstxy CX and at end after lastlocal layer',
 'maxiter': 400,
 'maxfun': 20000,
 'maxls': 20,
 'ftol': 1e-15,
 'gtol': 1e-10,
 'internal_seconds': 45,
 'external_seconds': 60,
 'case_seconds': 15,
 'optimizer_global_seconds': 35,
 'thread_count': 1,
 'numeric_tolerance': 1e-09,
 'restarts': 0,
 'ranking': 'deterministic reduced minimum-loss summary after both fixed-budget cases; no adaptive '
            'allocation',
 'adaptive_ranking': False,
 'early_stop_loss': 1e-22,
 'evidence_scope': 'Bounded reduced numerical diagnostics and fullV4 comparators; reduced success does '
                   'not imply fullV4 success; no exact certificate, incumbent promotion, topology '
                   'exclusion, stationary-point proof, adaptive allocation/deepening, or ten-CX lower '
                   'bound.'}
REQUIRED = set(FIXED) | {'script_sha256','family_sha256','checker_sha256','runtime_versions','dependency_sha256','plan_sha256','dependencies_file_sha256','input_manifest_sha256','frozen_inputs_sha256'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def versions():
    return {'python': '.'.join(map(str, sys.version_info[:3])), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'torch': torch.__version__, 'threadpoolctl': threadpoolctl.__version__}


def dependencies():
    return {'python': sha(sys.executable), 'numpy_init': sha(np.__file__),
            'numpy_multiarray': sha(np._core._multiarray_umath.__file__),
            'scipy_init': sha(scipy.__file__), 'scipy_lbfgsb': sha(scipy.optimize._lbfgsb.__file__),
            'torch_init': sha(torch.__file__), 'torch_C': sha(torch._C.__file__),
            'torch_cpu': sha(Path(torch.__file__).parent/'lib/libtorch_cpu.so'),
            'threadpoolctl': sha(threadpoolctl.__file__)}


def guard():
    if not HERE.name.isdigit() or HERE.parent.name!='runs' or HERE.parent.parent.name!='experiments':
        raise RuntimeError('reserved experiments/runs/<ID> workspace required')
    root,run_id=HERE.parents[2],int(HERE.name);raw=(HERE/'config.json').read_bytes();cfg=json.loads(raw)
    if raw!=json.dumps(cfg,sort_keys=True,separators=(',',':'),allow_nan=False).encode() or set(cfg)!=REQUIRED:
        raise RuntimeError('exact canonical raw config schema required')
    for k,v in FIXED.items():
        if cfg[k]!=v or type(cfg[k]) is not type(v):raise RuntimeError('fixed config mismatch '+k)
    for name,key in [('fit.py','script_sha256'),('family.py','family_sha256'),('PLAN.md','plan_sha256'),('dependencies.json','dependencies_file_sha256'),('input_manifest.json','input_manifest_sha256')]:
        if sha(HERE/name)!=cfg[key]:raise RuntimeError('frozen local hash mismatch '+name)
    if set(cfg['frozen_inputs_sha256'])!={'old_source.json','prefix.json','topologies.json'}:raise RuntimeError('local input set mismatch')
    for name,h in cfg['frozen_inputs_sha256'].items():
        if sha(HERE/name)!=h:raise RuntimeError('frozen local input mismatch '+name)
    if cfg['checker_sha256']!=sha(root/'check_circuit.py') or cfg['runtime_versions']!=versions() or cfg['dependency_sha256']!=dependencies():raise RuntimeError('checker/runtime/dependency mismatch')
    if json.loads((HERE/'dependencies.json').read_text())!={'runtime_versions':cfg['runtime_versions'],'dependency_sha256':cfg['dependency_sha256']}:raise RuntimeError('dependency snapshot mismatch')
    if (root/'code_snapshots'/cfg['script_sha256']).read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('reserved code snapshot mismatch')
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for rel,h in manifest['files'].items():
        if sha(root/rel)!=h:raise RuntimeError('immutable origin changed '+rel)
    source=json.loads((HERE/'old_source.json').read_text());prefix=json.loads((HERE/'prefix.json').read_text());plans=json.loads((HERE/'topologies.json').read_text())
    if sha(HERE/'old_source.json')!=cfg['source_sha256'] or source['n']!=3 or len(source['gates'])!=15 or [[g['control'],g['target']] for g in source['gates'] if g['gate']=='cx']!=cfg['old_topology'] or any(g['qubit'] not in (0,2) for g in source['gates'] if g['gate']!='cx'):raise RuntimeError('phase-correct old source convention mismatch')
    accepted=json.loads((root/'experiments/runs/411/candidate.json').read_text())
    if prefix['n']!=4 or prefix['gates']!=accepted['gates'][:6] or sum(g['gate']=='cx' for g in prefix['gates'])!=4:raise RuntimeError('accepted exact Bell/router prefix mismatch')
    expected_plans={'n':3,'wire_order':['x','u','y'],'old_topology':cfg['old_topology'],'new_topology':cfg['new_topology'],'rule':cfg['topology_rule'],'full_prefix_cx_count':4,'full_bare_cx_edges':[[3,2],[3,2]],'full_bare_cx_placement':['after first new xy CX, before local layer1','at end after local layer4'],'full_cx_count':10}
    if plans!=expected_plans:raise RuntimeError('complete topology/comparator rule differs')
    with sqlite3.connect(f"file:{root/'experiments.sqlite3'}?mode=ro",uri=True) as db:
        row=db.execute('SELECT status,code_sha256,config_json,method,agent,parent_id,target FROM attempts WHERE id=?',(run_id,)).fetchone()
    if row!=('running',cfg['script_sha256'],raw.decode(),cfg['method'],cfg['agent'],cfg['parent_run'],cfg['target']):raise RuntimeError('active exact raw-config ledger provenance mismatch')
    with sqlite3.connect(f"file:{root/'collaboration/board.sqlite3'}?mode=ro",uri=True) as db:
        board=db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?',(cfg['board_id'],)).fetchone()
    if board!=('running',cfg['agent'],run_id):raise RuntimeError('same-owner active board/run link mismatch')
    if any(os.environ[k]!='1' for k in THREAD_KEYS):raise RuntimeError('one-thread environment mismatch')
    return root,run_id,cfg,source,prefix,plans


def reduced_metrics(vectors,topology,family,tolerance):
    U=family.reduced_numpy(vectors,topology);A=U@family.target_v0_numpy()@U.conj().T;off=A-np.diag(np.diag(A))
    loss=float(np.vdot(off,off).real/8);unitarity=float(np.max(np.abs(U@U.conj().T-np.eye(8))))
    maxoff=float(np.max(np.abs(off)));root=float(np.max(np.abs(np.diag(A)**4-1)))
    return {'target':'Wv0=CZxy CXu-to-x; not V3','loss':loss,'unitarity_error':unitarity,'off_diagonal_error':maxoff,'eigenvalue_root_error':root,'valid_reduced_diagonalizer':max(unitarity,maxoff,root)<=tolerance,'reduced_cnot_count':4}


class BudgetStop(RuntimeError):
    pass


class TotalDeadline(BaseException):
    """Bypass case-level exception handling; preserve bests and stop the run."""
    pass


def main():
    started = time.monotonic()
    root,run_id,cfg,source,prefix,plans=guard()
    with (HERE/'execution_started.json').open('x') as marker:
        json.dump({'run_id': run_id, 'config_sha256': sha(HERE/'config.json'), 'started_unix': time.time()}, marker)
    remaining = cfg['internal_seconds']-(time.monotonic()-started)
    if remaining <= 0:
        raise BudgetStop('preflight exhausted total budget')
    def alarm(signum, frame):
        raise TotalDeadline('internal45-second deadline reached')
    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, remaining)
    states = []
    terminal_error = None
    try:
        sys.path.insert(0, str(root))
        import check_circuit as checker
        if Path(checker.__file__).resolve() != (root/'check_circuit.py').resolve():
            raise RuntimeError('loaded unexpected shared checker')
        spec = importlib.util.spec_from_file_location('reduced_v0_family', HERE/'family.py')
        family = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(family)
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        rng = np.random.Generator(np.random.PCG64(cfg['seed']))
        noises = [rng.normal(0.0, sigma, size=cfg['free_parameter_count']) for sigma in cfg['noise_sigmas']]
        # Allocate every seed draw at the start; early termination cannot change allocations.
        write(HERE/'all_noise30.json', [noise.tolist() for noise in noises])
        cases_dir = HERE/'cases'
        cases_dir.mkdir(exist_ok=False)
        with threadpool_limits(limits=cfg['thread_count']):
            base,old_topology,matrices=family.collapse(source,checker)
            if old_topology!=cfg['old_topology'] or np.shape(base)!=(30,):raise RuntimeError('source collapse shape/topology differs')
            write(HERE/'base30.json',{'base30':base.tolist(),'old_topology':old_topology,'new_topology':cfg['new_topology']})
            write(HERE/'old_base_candidate.json',family.serialize_reduced(base,old_topology))
            write(HERE/'old_local_matrices.json',[[[[[float(z.real),float(z.imag)] for z in row] for row in U] for U in layer] for layer in matrices])
            # Retain both fixed perturbation points before either optimizer starts.
            for ordinal,(sigma,noise) in enumerate(zip(cfg['noise_sigmas'],noises)):
                path=cases_dir/f'{ordinal:02d}_sigma{sigma:g}';path.mkdir();topology=cfg['new_topology'];initial=base+noise
                for label,point in [('base',base),('initial',initial),('best',initial)]:
                    reduced=family.serialize_reduced(point,topology);full=family.serialize_full(point,topology,prefix)
                    write(path/(label+'_reduced_candidate.json'),reduced);write(path/(label+'_full_candidate.json'),full)
                write(path/'vectors.json',{'base30':base.tolist(),'noise30':noise.tolist(),'initial30':initial.tolist(),'old_topology':old_topology,'new_topology':topology,'noise_sigma':sigma})
                base_metric=reduced_metrics(base,topology,family,cfg['numeric_tolerance']);initial_metric=reduced_metrics(initial,topology,family,cfg['numeric_tolerance'])
                base_full=checker.evaluate(family.serialize_full(base,topology,prefix),tolerance=cfg['numeric_tolerance']);initial_full=checker.evaluate(family.serialize_full(initial,topology,prefix),tolerance=cfg['numeric_tolerance'])
                if base_full['cnot_count']!=10 or initial_full['cnot_count']!=10:raise RuntimeError('full comparator count mismatch')
                write(path/'base_metrics.json',{'reduced':base_metric,'full_v4':base_full});write(path/'initial_metrics.json',{'reduced':initial_metric,'full_v4':initial_full})
                states.append({'ordinal':ordinal,'path':path,'noise_sigma':sigma,'topology':topology,'initial':initial,'base_reduced':base_metric,'initial_reduced':initial_metric,'base_full_checker':base_full,'initial_full_checker':initial_full,'initial_loss':initial_metric['loss'],'best_x':initial.copy(),'best_loss':initial_metric['loss'],'best_gradient':None,'best_evaluation':None,'calls':0,'iterations':0,'history':[],'termination':'not_started','optimizer':None,'case_seconds':0.})
            global_deadline = started+cfg['optimizer_global_seconds']
            for state in states:
                case_started = time.monotonic()
                if case_started >= global_deadline:
                    state['termination'] = 'skipped_global_budget'
                    continue
                deadline = min(case_started+cfg['case_seconds'], global_deadline)
                word = family.TorchWord(state['topology'])
                def objective(z):
                    if time.monotonic() >= deadline:
                        raise BudgetStop('case15-second or global optimization deadline')
                    if state['calls'] >= cfg['maxfun']:
                        raise BudgetStop('max20000 objective/gradient calls')
                    x = np.asarray(z, dtype=np.float64).copy()
                    variable = torch.tensor(x, dtype=torch.float64, requires_grad=True)
                    loss_tensor = word.loss(variable)
                    grad = torch.autograd.grad(loss_tensor, variable)[0].detach().numpy()
                    loss = float(loss_tensor.detach())
                    if not np.isfinite(loss) or not np.isfinite(grad).all():
                        raise FloatingPointError('nonfinite objective/gradient')
                    state['calls'] += 1
                    state['history'].append({'evaluation': state['calls'], 'loss': loss, 'gradient_norm': float(np.linalg.norm(grad))})
                    if loss <= state['best_loss']:
                        state.update(best_x=x, best_loss=loss, best_gradient=grad.copy(), best_evaluation=state['calls'])
                    if loss <= cfg['early_stop_loss']:
                        raise BudgetStop('numerical loss trigger1e-22; exact certificate remains separate')
                    return loss, grad
                def callback(x):
                    state['iterations'] += 1
                    write(state['path']/'best_reduced_candidate.json',family.serialize_reduced(state['best_x'],state['topology']))
                    write(state['path']/'best_full_candidate.json',family.serialize_full(state['best_x'],state['topology'],prefix))
                    write(state['path']/'best_vector_checkpoint.json', {'best30': state['best_x'].tolist(), 'loss': state['best_loss'], 'evaluation': state['best_evaluation'], 'iterations': state['iterations']})
                try:
                    r = minimize(objective, state['initial'], jac=True, method='L-BFGS-B', callback=callback,
                                 options={key: cfg[key] for key in ('maxiter', 'maxfun', 'maxls', 'ftol', 'gtol')})
                    state['termination'] = 'optimizer_returned: '+str(r.message)
                    state['optimizer'] = {'success': bool(r.success), 'message': str(r.message), 'nit': int(r.nit), 'nfev': int(r.nfev), 'njev': int(r.njev), 'fun': float(r.fun)}
                    state['iterations'] = int(r.nit)
                except BudgetStop as error:
                    state['termination'] = 'budget_or_loss_stop: '+str(error)
                except Exception as error:
                    state['termination'] = f'case_exception: {type(error).__name__}: {error}'
                state['case_seconds'] = time.monotonic()-case_started
                write(state['path']/'best_reduced_candidate.json',family.serialize_reduced(state['best_x'],state['topology']))
                write(state['path']/'best_full_candidate.json',family.serialize_full(state['best_x'],state['topology'],prefix))
            records=[]
            for state in states:
                reduced=family.serialize_reduced(state['best_x'],state['topology']);full=family.serialize_full(state['best_x'],state['topology'],prefix)
                write(state['path']/'best_reduced_candidate.json',reduced);write(state['path']/'best_full_candidate.json',full)
                reduced_checked=reduced_metrics(state['best_x'],state['topology'],family,cfg['numeric_tolerance']);full_checked=checker.evaluate(full,tolerance=cfg['numeric_tolerance'])
                if len(reduced['gates'])!=14 or len(full['gates'])!=22 or full_checked['cnot_count']!=10:raise RuntimeError('reduced/full endpoint count mismatch')
                reduced_file=str((state['path']/'best_reduced_candidate.json').relative_to(HERE));full_file=str((state['path']/'best_full_candidate.json').relative_to(HERE))
                record={'case_index':state['ordinal'],'noise_sigma':state['noise_sigma'],'old_topology':plans['old_topology'],'topology':state['topology'],'free_parameter_count':30,'base_reduced':state['base_reduced'],'initial_reduced':state['initial_reduced'],'best_reduced':reduced_checked,'base_full_checker':state['base_full_checker'],'initial_full_checker':state['initial_full_checker'],'best_full_checker':full_checked,'initial_loss':state['initial_loss'],'best_loss':state['best_loss'],'best30':state['best_x'].tolist(),'best_gradient':None if state['best_gradient'] is None else state['best_gradient'].tolist(),'best_evaluation':state['best_evaluation'],'evaluations':state['calls'],'iterations':state['iterations'],'history':state['history'],'termination':state['termination'],'optimizer':state['optimizer'],'optimizer_seconds':state['case_seconds'],'reduced_candidate_file':reduced_file,'full_candidate_file':full_file,'reduced_candidate_sha256':sha(HERE/reduced_file),'full_candidate_sha256':sha(HERE/full_file)}
                write(state['path']/'result.json',record);records.append(record)
            best=min(records,key=lambda r:(r['best_loss'],r['case_index']))
            write(HERE/'best_reduced_candidate.json',json.loads((HERE/best['reduced_candidate_file']).read_text()))
            write(HERE/'best_full_candidate.json',json.loads((HERE/best['full_candidate_file']).read_text()))
            result={'status':'completed','run_id':run_id,'source_run':411,'source_sha256':cfg['source_sha256'],'script_sha256':cfg['script_sha256'],'family_sha256':cfg['family_sha256'],'config_sha256':sha(HERE/'config.json'),'input_manifest_sha256':cfg['input_manifest_sha256'],'runtime_versions':versions(),'dependency_sha256':cfg['dependency_sha256'],'thread_environment':{k:os.environ[k] for k in THREAD_KEYS},'torch_threads':torch.get_num_threads(),'torch_interop_threads':torch.get_num_interop_threads(),'noise_sha256':sha(HERE/'all_noise30.json'),'case_count':len(records),'records':records,'old_topology':plans['old_topology'],'new_topology':plans['new_topology'],'best':best,'best_reduced_candidate_sha256':sha(HERE/'best_reduced_candidate.json'),'best_full_candidate_sha256':sha(HERE/'best_full_candidate.json'),'seconds':time.monotonic()-started,'ranking':cfg['ranking'],'adaptive_ranking':False,'evidence_scope':cfg['evidence_scope']}
            write(HERE/'result.json',result)
            print(json.dumps({'run_id':run_id,'case_count':len(records),'best_case_index':best['case_index'],'records':[{'case_index':r['case_index'],'reduced_loss':r['best_loss'],'reduced_maxoff':r['best_reduced']['off_diagonal_error'],'reduced_pass':r['best_reduced']['valid_reduced_diagonalizer'],'full_maxoff':r['best_full_checker']['off_diagonal_error'],'full_pass':r['best_full_checker']['valid_diagonalizer']} for r in records],'seconds':result['seconds']},sort_keys=True))
    except BaseException as error:
        terminal_error = {'type': type(error).__name__, 'message': str(error)}
        # Preserve best vectors and serialize gate lists; no objective/checker call.
        signal.setitimer(signal.ITIMER_REAL, 0)
        for state in states:
            write(state['path']/'best_reduced_candidate.json',family.serialize_reduced(state['best_x'],state['topology']))
            write(state['path']/'best_full_candidate.json',family.serialize_full(state['best_x'],state['topology'],prefix))
            write(state['path']/'terminal_state.json', {'best30': state['best_x'].tolist(), 'best_loss': state['best_loss'],
                  'calls': state['calls'], 'iterations': state['iterations'], 'history': state['history'],
                  'termination': state['termination'], 'terminal_error': terminal_error})
        write(HERE/'failure.json', {'run_id': run_id, 'terminal_error': terminal_error, 'prepared_cases': len(states), 'seconds': time.monotonic()-started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
