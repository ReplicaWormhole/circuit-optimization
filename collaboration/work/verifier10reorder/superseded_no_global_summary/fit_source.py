#!/usr/bin/env python3
"""Two frozen reordered-topology diagnostics, requiring coordinator reservation."""
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
 'purpose': 'two bounded full-local ten-CX reordered-topology diagnostics from immutable run416 '
            'deletion7 best vector',
 'method': 'SciPy L-BFGS-B with Torch analytic gradients; full132 SU2 locals; one reordered ten-CX '
           'topology and two frozen perturbations',
 'board_id': 148,
 'agent': 'numerical10reorder',
 'parent_run': 416,
 'target': 'V4',
 'source_filename': 'source_case_result.json',
 'topology_filename': 'topologies.json',
 'source_case_result_sha256': '853d53241573363a5baf9d86f833cdf9ec5c6304cc0b2b2fa168466015864ef7',
 'source_best_candidate_sha256': '26b9f34ab733deb77127dcd666929ae8c8a515746892f5ddced57fdf4d36985b',
 'seed': 10012301,
 'rng': 'NumPy Generator PCG64',
 'noise_sigmas': [0.025, 0.25],
 'noise_rule': 'one PCG64(seed) generator; two consecutive normal(0,sigma,size=132) draws in case '
               'order0..1 with sigma0.025 then0.25; allocate all upfront',
 'case_count': 2,
 'free_parameter_count': 132,
 'local_layer_count': 11,
 'expected_cx_count': 10,
 'topology_rule': 'swap old CX ordinals6 and7; retain all other edges and all layer coordinates; this '
                  'changes noncommuting CX order',
 'new_topology': [[0, 2], [1, 3], [0, 1], [2, 3], [0, 2], [3, 2], [1, 0], [0, 2], [1, 2], [3, 2]],
 'old_topology': [[0, 2], [1, 3], [0, 1], [2, 3], [0, 2], [3, 2], [0, 2], [1, 0], [1, 2], [3, 2]],
 'maxiter': 250,
 'maxfun': 10000,
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
 'ranking': False,
 'early_stop_loss': 1e-22,
 'evidence_scope': 'Two bounded numerical diagnostics; no exact certificate, incumbent promotion, '
                   'topology exclusion, stationary-point proof, ranking of starts, or ten-CX lower '
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
    if not HERE.name.isdigit() or HERE.parent.name != 'runs' or HERE.parent.parent.name != 'experiments':
        raise RuntimeError('execute only from experiments/runs/<reserved-id>')
    root, run_id = HERE.parents[2], int(HERE.name)
    config_raw = (HERE/'config.json').read_bytes()
    cfg = json.loads(config_raw)
    if config_raw != json.dumps(cfg, sort_keys=True, separators=(',', ':'), allow_nan=False).encode():
        raise RuntimeError('config must be exact canonical raw JSON bytes')
    if set(cfg) != REQUIRED:
        raise RuntimeError('config-key mismatch')
    for key, expected in FIXED.items():
        if type(cfg[key]) is not type(expected) or cfg[key] != expected:
            raise RuntimeError('frozen config value mismatch: '+key)
    for name,key in [('fit.py','script_sha256'),('family.py','family_sha256'),('PLAN.md','plan_sha256'),('dependencies.json','dependencies_file_sha256'),('input_manifest.json','input_manifest_sha256')]:
        if sha(HERE/name) != cfg[key]:
            raise RuntimeError('frozen local hash mismatch: '+name)
    expected_inputs = {'source_case_result.json','source_best_candidate.json','topologies.json','base_vector.json','base_candidate.json'}
    if set(cfg['frozen_inputs_sha256']) != expected_inputs:
        raise RuntimeError('frozen local input set mismatch')
    for name,expected in cfg['frozen_inputs_sha256'].items():
        if sha(HERE/name) != expected:
            raise RuntimeError('frozen local input hash mismatch: '+name)
    if cfg['checker_sha256'] != sha(root/'check_circuit.py') or cfg['runtime_versions'] != versions() or cfg['dependency_sha256'] != dependencies():
        raise RuntimeError('checker/runtime/dependency mismatch')
    if json.loads((HERE/'dependencies.json').read_text()) != {'runtime_versions': cfg['runtime_versions'], 'dependency_sha256': cfg['dependency_sha256']}:
        raise RuntimeError('dependency snapshot content mismatch')
    if (root/'code_snapshots'/cfg['script_sha256']).read_bytes() != Path(__file__).read_bytes():
        raise RuntimeError('reserved code snapshot mismatch')
    manifest = json.loads((HERE/'input_manifest.json').read_text())
    for relative,expected in manifest['files'].items():
        if sha(root/relative) != expected:
            raise RuntimeError('immutable source input changed: '+relative)
    source = json.loads((HERE/'source_case_result.json').read_text())
    old_word = json.loads((HERE/'source_best_candidate.json').read_text())
    if sha(HERE/'source_case_result.json') != cfg['source_case_result_sha256'] or sha(HERE/'source_best_candidate.json') != cfg['source_best_candidate_sha256']:
        raise RuntimeError('immutable copied source hash mismatch')
    if source['delete_cx_ordinal'] != 7 or source['topology'] != cfg['old_topology'] or len(source['best132']) != 132 or not all(math.isfinite(x) for x in source['best132']):
        raise RuntimeError('immutable source coordinates/topology mismatch')
    base = json.loads((HERE/'base_vector.json').read_text())
    if base != {'base132': source['best132'], 'source_run':416, 'source_deletion_ordinal':7, 'source_case_result_sha256':cfg['source_case_result_sha256']}:
        raise RuntimeError('base coordinates were not transplanted unchanged')
    plans = json.loads((HERE/'topologies.json').read_text())
    if plans != {'old_topology':cfg['old_topology'],'new_topology':cfg['new_topology'],'rule':cfg['topology_rule'],'old_to_new_cx_ordinal':[0,1,2,3,4,5,7,6,8,9]}:
        raise RuntimeError('complete topology rule differs')
    literal = json.loads((HERE/'base_candidate.json').read_text())
    expected_gates = [dict(gate) for gate in old_word['gates']]
    index = 0
    for gate in expected_gates:
        if gate['gate'] == 'cx':
            gate['control'],gate['target'] = cfg['new_topology'][index]
            index += 1
    if literal['n'] != 4 or literal['gates'] != expected_gates or index != 10 or len(expected_gates) != 54:
        raise RuntimeError('literal new base gate list differs')
    with sqlite3.connect(f"file:{root/'experiments.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute('SELECT status,code_sha256,config_json,method,agent,parent_id,target FROM attempts WHERE id=?',(run_id,)).fetchone()
        parent = db.execute('SELECT status,agent,parent_id,target FROM attempts WHERE id=416').fetchone()
    if row is None or row[0] != 'running' or row[1] != cfg['script_sha256'] or row[2].encode() != config_raw or row[3] != cfg['method'] or row[4:] != (cfg['agent'],cfg['parent_run'],cfg['target']):
        raise RuntimeError('active ledger code/rawconfig/method/agent/parent/target mismatch')
    if parent != ('failed','numerical10',411,'V4'):
        raise RuntimeError('immutable source run status/provenance mismatch')
    with sqlite3.connect(f"file:{root/'collaboration/board.sqlite3'}?mode=ro",uri=True) as db:
        board = db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?',(cfg['board_id'],)).fetchone()
    if board != ('running',cfg['agent'],run_id):
        raise RuntimeError('active same-owner board/run link mismatch')
    if any(os.environ[key] != '1' for key in THREAD_KEYS):
        raise RuntimeError('thread environment mismatch')
    return root,run_id,cfg,source,plans


class BudgetStop(RuntimeError):
    pass


class TotalDeadline(BaseException):
    """Bypass case-level exception handling; preserve bests and stop the run."""
    pass


def main():
    started = time.monotonic()
    root, run_id, cfg, source, plans = guard()
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
        spec = importlib.util.spec_from_file_location('ten_full_local_family', HERE/'family.py')
        family = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(family)
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        rng = np.random.Generator(np.random.PCG64(cfg['seed']))
        noises = [rng.normal(0.0, sigma, size=cfg['free_parameter_count']) for sigma in cfg['noise_sigmas']]
        # Allocate every seed draw at the start; early termination cannot change allocations.
        write(HERE/'all_noise132.json', [noise.tolist() for noise in noises])
        cases_dir = HERE/'cases'
        cases_dir.mkdir(exist_ok=False)
        with threadpool_limits(limits=cfg['thread_count']):
            # Save all complete initial words and vectors before any optimizer starts.
            for ordinal, (sigma,noise) in enumerate(zip(cfg['noise_sigmas'],noises)):
                path = cases_dir/f'{ordinal:02d}_sigma{sigma:g}'
                path.mkdir()
                base = np.asarray(source['best132'],dtype=np.float64)
                topology = plans['new_topology']
                initial = base+noise
                base_word, initial_word = family.serialize(base,topology),family.serialize(initial,topology)
                for name,word in [('base_candidate.json',base_word),('literal_base_candidate.json',json.loads((HERE/'base_candidate.json').read_text())),('initial_candidate.json',initial_word),('best_candidate.json',initial_word)]:
                    write(path/name,word)
                write(path/'vectors.json',{'base132':base.tolist(),'noise132':noise.tolist(),'initial132':initial.tolist(),'old_topology':plans['old_topology'],'new_topology':topology,'noise_sigma':sigma})
                base_checked = checker.evaluate(base_word,tolerance=cfg['numeric_tolerance'])
                base_U = family.unitary_numpy(base,topology)
                base_A = base_U @ family.right_shift_numpy() @ base_U.conj().T
                base_off = base_A-np.diag(np.diag(base_A))
                base_loss = float(np.vdot(base_off,base_off).real/16)
                write(path/'base_metrics.json',{'loss':base_loss,'checker':base_checked})
                initial_checked = checker.evaluate(initial_word, tolerance=cfg['numeric_tolerance'])
                initial_U = family.unitary_numpy(initial, topology)
                A = initial_U @ family.right_shift_numpy() @ initial_U.conj().T
                off = A-np.diag(np.diag(A))
                initial_loss = float(np.vdot(off, off).real/16)
                if initial_checked['cnot_count'] != cfg['expected_cx_count']:
                    raise RuntimeError('serialized initial CX count differs from ten')
                state = {'ordinal': ordinal, 'path': path, 'noise_sigma':sigma, 'topology': topology,
                         'base_loss':base_loss,'base_checker':base_checked,'initial': initial, 'initial_checker': initial_checked, 'initial_loss': initial_loss,
                         'best_x': initial.copy(), 'best_loss': initial_loss, 'best_gradient': None,
                         'best_evaluation': None, 'calls': 0, 'iterations': 0, 'history': [],
                         'termination': 'not_started', 'optimizer': None, 'case_seconds': 0.0}
                states.append(state)
                write(path/'initial_metrics.json', {'loss': initial_loss, 'checker': initial_checked})
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
                        raise BudgetStop('max10000 objective/gradient calls')
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
                    write(state['path']/'best_candidate.json', family.serialize(state['best_x'], state['topology']))
                    write(state['path']/'best_vector_checkpoint.json', {'best132': state['best_x'].tolist(), 'loss': state['best_loss'], 'evaluation': state['best_evaluation'], 'iterations': state['iterations']})
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
                write(state['path']/'best_candidate.json', family.serialize(state['best_x'], state['topology']))
            records = []
            for state in states:
                best_word = family.serialize(state['best_x'], state['topology'])
                write(state['path']/'best_candidate.json', best_word)
                checked = checker.evaluate(best_word, tolerance=cfg['numeric_tolerance'])
                if checked['cnot_count'] != cfg['expected_cx_count']:
                    raise RuntimeError('best candidate CX count differs from ten')
                record = {'case_index':state['ordinal'],'noise_sigma':state['noise_sigma'],'old_topology':plans['old_topology'],'topology':state['topology'],
                          'free_parameter_count':cfg['free_parameter_count'],'base_loss':state['base_loss'],'base_checker':state['base_checker'],'initial_loss':state['initial_loss'],
                          'initial_checker': state['initial_checker'], 'best_loss': state['best_loss'], 'best_checker': checked,
                          'best132': state['best_x'].tolist(), 'best_gradient': None if state['best_gradient'] is None else state['best_gradient'].tolist(),
                          'best_evaluation': state['best_evaluation'], 'evaluations': state['calls'], 'iterations': state['iterations'],
                          'history': state['history'], 'termination': state['termination'], 'optimizer': state['optimizer'],
                          'optimizer_seconds': state['case_seconds'], 'candidate_file': str((state['path']/'best_candidate.json').relative_to(HERE)),
                          'candidate_sha256': sha(state['path']/'best_candidate.json')}
                write(state['path']/'result.json', record)
                records.append(record)
            result = {'status':'completed','run_id':run_id,'source_run':416,'source_case_result_sha256':cfg['source_case_result_sha256'],
                      'script_sha256':cfg['script_sha256'],'family_sha256':cfg['family_sha256'],'config_sha256':sha(HERE/'config.json'),
                      'input_manifest_sha256':cfg['input_manifest_sha256'],'runtime_versions':versions(),'dependency_sha256':cfg['dependency_sha256'],
                      'thread_environment':{key:os.environ[key] for key in THREAD_KEYS},'torch_threads':torch.get_num_threads(),
                      'torch_interop_threads':torch.get_num_interop_threads(),'noise_sha256':sha(HERE/'all_noise132.json'),
                      'case_count':len(records),'records':records,'old_topology':plans['old_topology'],'new_topology':plans['new_topology'],
                      'seconds':time.monotonic()-started,'ranking':False,'evidence_scope':cfg['evidence_scope']}
            write(HERE/'result.json',result)
            print(json.dumps({'run_id':run_id,'case_count':len(records),'records':[{'case_index':r['case_index'],'noise_sigma':r['noise_sigma'],'best_loss':r['best_loss'],'maxoff':r['best_checker']['off_diagonal_error'],'pass':r['best_checker']['valid_diagonalizer']} for r in records],'seconds':result['seconds']},sort_keys=True))
    except BaseException as error:
        terminal_error = {'type': type(error).__name__, 'message': str(error)}
        # Preserve best vectors and serialize gate lists; no objective/checker call.
        signal.setitimer(signal.ITIMER_REAL, 0)
        for state in states:
            write(state['path']/'best_candidate.json', family.serialize(state['best_x'], state['topology']))
            write(state['path']/'terminal_state.json', {'best132': state['best_x'].tolist(), 'best_loss': state['best_loss'],
                  'calls': state['calls'], 'iterations': state['iterations'], 'history': state['history'],
                  'termination': state['termination'], 'terminal_error': terminal_error})
        write(HERE/'failure.json', {'run_id': run_id, 'terminal_error': terminal_error, 'prepared_cases': len(states), 'seconds': time.monotonic()-started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
