#!/usr/bin/env python3
"""One frozen twelve-case deletion refit, requiring coordinator reservation."""
from __future__ import annotations
import os
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
for key in THREAD_KEYS:
    os.environ[key] = '1'
import hashlib
import importlib.util
import json
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
FIXED = {
    'schema_version': 1, 'purpose': 'one bounded full-local single-deletion refit of the certified run406 exact12 circuit',
    'method': 'SciPy L-BFGS-B with Torch analytic gradients; full144 SU2 locals; twelve exact12 single-CX deletions',
    'board_id': 132, 'source_filename': 'source_candidate.json', 'topology_filename': 'topologies.json',
    'source_candidate_sha256': '58a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2',
    'seed': 10012101, 'rng': 'NumPy Generator PCG64', 'noise_sigma': 0.025,
    'noise_rule': 'one PCG64(seed) generator; twelve consecutive normal(0,sigma,size=144) draws in deletion ordinal order0..11',
    'case_count': 12, 'free_parameter_count': 144, 'local_layer_count': 12, 'expected_cx_count': 11,
    'topology_rule': 'delete each chronological CX ordinal0..11 once; retain all other CX order/direction; twelve free4wire SU2 layers',
    'maxiter': 80, 'maxfun': 4000, 'maxls': 20, 'ftol': 1e-15, 'gtol': 1e-10,
    'internal_seconds': 110, 'external_seconds': 120, 'case_seconds': 8,
    'optimizer_global_seconds': 102, 'thread_count': 1, 'numeric_tolerance': 1e-9,
    'restarts': 0, 'early_stop_loss': 1e-22,
    'evidence_scope': 'One bounded numerical diagnostic; no exact certificate, incumbent promotion, topology exclusion, stationary-point proof, or eleven-CX lower bound.',
}
REQUIRED = set(FIXED) | {'script_sha256', 'family_sha256', 'topologies_sha256', 'checker_sha256', 'runtime_versions', 'dependency_sha256'}


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


def topology_records(candidate):
    indexed = [(position, gate['control'], gate['target']) for position, gate in enumerate(candidate['gates']) if gate['gate'] == 'cx']
    if candidate['n'] != 4 or len(indexed) != 12:
        raise ValueError('source must be the complete four-qubit twelve-CX word')
    return [{'delete_cx_ordinal': index, 'source_gate_position': position, 'deleted_edge': [control, target],
             'topology': [[c, t] for j, (p, c, t) in enumerate(indexed) if j != index]}
            for index, (position, control, target) in enumerate(indexed)]


def guard():
    if not HERE.name.isdigit() or HERE.parent.name != 'runs' or HERE.parent.parent.name != 'experiments':
        raise RuntimeError('execute only from experiments/runs/<reserved-id>')
    root, run_id = HERE.parents[2], int(HERE.name)
    cfg = json.loads((HERE/'config.json').read_text())
    if set(cfg) != REQUIRED:
        raise RuntimeError(f'config-key mismatch: {sorted(set(cfg)^REQUIRED)}')
    for key, expected in FIXED.items():
        if type(cfg[key]) is not type(expected) or cfg[key] != expected:
            raise RuntimeError(f'frozen config value mismatch: {key}')
    if cfg['script_sha256'] != sha(__file__) or cfg['family_sha256'] != sha(HERE/'family.py'):
        raise RuntimeError('script/family hash mismatch')
    if cfg['checker_sha256'] != sha(root/'check_circuit.py'):
        raise RuntimeError('checker source hash mismatch')
    if cfg['source_candidate_sha256'] != sha(HERE/cfg['source_filename']) or cfg['topologies_sha256'] != sha(HERE/cfg['topology_filename']):
        raise RuntimeError('source candidate/topology hash mismatch')
    if cfg['runtime_versions'] != versions() or cfg['dependency_sha256'] != dependencies():
        raise RuntimeError('runtime version/dependency hash mismatch')
    candidate = json.loads((HERE/cfg['source_filename']).read_text())
    plans = json.loads((HERE/cfg['topology_filename']).read_text())
    if plans != topology_records(candidate):
        raise RuntimeError('topology serialization differs from the exact single-deletion rule')
    with sqlite3.connect(f"file:{root/'experiments.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute('SELECT status,code_sha256,config_json,method FROM attempts WHERE id=?', (run_id,)).fetchone()
    if row is None or row[0] != 'running' or row[1] != cfg['script_sha256'] or json.loads(row[2]) != cfg or row[3] != cfg['method']:
        raise RuntimeError('active ledger code/config/method mismatch')
    with sqlite3.connect(f"file:{root/'collaboration/board.sqlite3'}?mode=ro", uri=True) as db:
        board = db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?', (cfg['board_id'],)).fetchone()
    if board != ('running', 'numerical_sol', run_id):
        raise RuntimeError('active board ownership/run link mismatch')
    if any(os.environ[key] != '1' for key in THREAD_KEYS):
        raise RuntimeError('thread environment mismatch')
    return root, run_id, cfg, candidate, plans


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
        raise TotalDeadline('internal110-second deadline reached')
    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, remaining)
    states = []
    terminal_error = None
    try:
        sys.path.insert(0, str(root))
        import check_circuit as checker
        if Path(checker.__file__).resolve() != (root/'check_circuit.py').resolve():
            raise RuntimeError('loaded unexpected shared checker')
        spec = importlib.util.spec_from_file_location('eleven_full_local_family', HERE/'family.py')
        family = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(family)
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        rng = np.random.Generator(np.random.PCG64(cfg['seed']))
        noises = [rng.normal(0.0, cfg['noise_sigma'], size=cfg['free_parameter_count']) for case in plans]
        # Allocate every seed draw at the start; early termination cannot change allocations.
        write(HERE/'all_noise144.json', [noise.tolist() for noise in noises])
        cases_dir = HERE/'cases'
        cases_dir.mkdir(exist_ok=False)
        with threadpool_limits(limits=cfg['thread_count']):
            # Save all complete initial words and vectors before any optimizer starts.
            for plan, noise in zip(plans, noises):
                ordinal = plan['delete_cx_ordinal']
                path = cases_dir/f'{ordinal:02d}_delete_cx{ordinal:02d}'
                path.mkdir()
                base, topology, deleted, matrices, position = family.collapse(source, ordinal, checker)
                if topology != plan['topology'] or position != plan['source_gate_position']:
                    raise RuntimeError('collapse topology differs from frozen rule')
                initial = base+noise
                base_word, initial_word = family.serialize(base, topology), family.serialize(initial, topology)
                for name, word in [('deleted_source.json', deleted), ('base_candidate.json', base_word), ('initial_candidate.json', initial_word), ('best_candidate.json', initial_word)]:
                    write(path/name, word)
                write(path/'vectors.json', {'base144': base.tolist(), 'noise144': noise.tolist(), 'initial144': initial.tolist()})
                write(path/'base_local_matrices.json', [[[[[float(z.real), float(z.imag)] for z in row] for row in U] for U in layer] for layer in matrices])
                initial_checked = checker.evaluate(initial_word, tolerance=cfg['numeric_tolerance'])
                initial_U = family.unitary_numpy(initial, topology)
                A = initial_U @ family.right_shift_numpy() @ initial_U.conj().T
                off = A-np.diag(np.diag(A))
                initial_loss = float(np.vdot(off, off).real/16)
                if initial_checked['cnot_count'] != cfg['expected_cx_count']:
                    raise RuntimeError('serialized initial CX count differs from eleven')
                state = {'ordinal': ordinal, 'path': path, 'plan': plan, 'topology': topology,
                         'initial': initial, 'initial_checker': initial_checked, 'initial_loss': initial_loss,
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
                        raise BudgetStop('case8-second or global optimization deadline')
                    if state['calls'] >= cfg['maxfun']:
                        raise BudgetStop('max4000 objective/gradient calls')
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
                    write(state['path']/'best_vector_checkpoint.json', {'best144': state['best_x'].tolist(), 'loss': state['best_loss'], 'evaluation': state['best_evaluation'], 'iterations': state['iterations']})
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
                    raise RuntimeError('best candidate CX count differs from eleven')
                record = {'delete_cx_ordinal': state['ordinal'], 'source_gate_position': state['plan']['source_gate_position'],
                          'deleted_edge': state['plan']['deleted_edge'], 'topology': state['topology'],
                          'free_parameter_count': cfg['free_parameter_count'], 'initial_loss': state['initial_loss'],
                          'initial_checker': state['initial_checker'], 'best_loss': state['best_loss'], 'best_checker': checked,
                          'best144': state['best_x'].tolist(), 'best_gradient': None if state['best_gradient'] is None else state['best_gradient'].tolist(),
                          'best_evaluation': state['best_evaluation'], 'evaluations': state['calls'], 'iterations': state['iterations'],
                          'history': state['history'], 'termination': state['termination'], 'optimizer': state['optimizer'],
                          'optimizer_seconds': state['case_seconds'], 'candidate_file': str((state['path']/'best_candidate.json').relative_to(HERE)),
                          'candidate_sha256': sha(state['path']/'best_candidate.json')}
                write(state['path']/'result.json', record)
                records.append(record)
            best = min(records, key=lambda r: (r['best_loss'], r['delete_cx_ordinal']))
            write(HERE/'best_candidate.json', json.loads((HERE/best['candidate_file']).read_text()))
            result = {'status': 'completed', 'run_id': run_id, 'source_candidate_sha256': cfg['source_candidate_sha256'],
                      'script_sha256': cfg['script_sha256'], 'family_sha256': cfg['family_sha256'], 'config_sha256': sha(HERE/'config.json'),
                      'runtime_versions': versions(), 'dependency_sha256': cfg['dependency_sha256'],
                      'thread_environment': {key: os.environ[key] for key in THREAD_KEYS}, 'torch_threads': torch.get_num_threads(),
                      'torch_interop_threads': torch.get_num_interop_threads(), 'noise_sha256': sha(HERE/'all_noise144.json'),
                      'case_count': len(records), 'best': best, 'records': records,
                      'best_candidate_sha256': sha(HERE/'best_candidate.json'), 'seconds': time.monotonic()-started,
                      'evidence_scope': cfg['evidence_scope']}
            write(HERE/'result.json', result)
            print(json.dumps({'run_id': run_id, 'case_count': len(records), 'best_loss': best['best_loss'],
                              'best_maxoff': best['best_checker']['off_diagonal_error'], 'pass': best['best_checker']['valid_diagonalizer'],
                              'seconds': result['seconds']}, sort_keys=True))
    except BaseException as error:
        terminal_error = {'type': type(error).__name__, 'message': str(error)}
        # Checkpoint the current best vectors and words without another matrix call.
        signal.setitimer(signal.ITIMER_REAL, 0)
        for state in states:
            write(state['path']/'best_candidate.json', family.serialize(state['best_x'], state['topology']))
            write(state['path']/'terminal_state.json', {'best144': state['best_x'].tolist(), 'best_loss': state['best_loss'],
                  'calls': state['calls'], 'iterations': state['iterations'], 'history': state['history'],
                  'termination': state['termination'], 'terminal_error': terminal_error})
        write(HERE/'failure.json', {'run_id': run_id, 'terminal_error': terminal_error, 'prepared_cases': len(states), 'seconds': time.monotonic()-started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
