#!/usr/bin/env python3
"""One separately reserved SVD Newton refinement of run409 deletion9 only."""
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
import torch
import threadpoolctl
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
FIXED = {
    'schema_version': 1,
    'purpose': 'one fresh minimum-norm Newton precision refinement of run409 deletion9 eleven-CX near-hit',
    'method': 'SVD minimum-norm Newton on480 real-imag off-diagonal residuals with full144 free SU2 locals',
    'board_id': 136, 'parent_run': 409, 'deletion_ordinal': 9,
    'seed': 0, 'random_draws': 0, 'restarts': 0,
    'source_point_filename': 'source_point.json', 'initial_candidate_filename': 'initial_candidate.json',
    'expected_cx_count': 11, 'parameter_count': 144, 'residual_count': 480,
    'svd_rcond': 1e-8, 'max_steps': 20, 'max_line_search': 8,
    'alpha_rule': 'alpha=2**(-j), j=0..7; accept strict decrease of squared residual norm',
    'fd_step': 1e-6, 'fd_absolute_tolerance': 1e-6, 'fd_relative_tolerance': 1e-5,
    'fd_direction_rule': 'cos(1),cos(2),...,cos(144), normalized by Euclidean norm; no RNG',
    'target_maxoff': 1e-13, 'checker_tolerance': 1e-9,
    'internal_seconds': 55, 'external_seconds': 60, 'thread_count': 1,
    'evidence_scope': 'Numerical precision refinement only; exact recognition/certification is separate. No optimum, local minimum, topology exclusion, or lower-bound claim.',
}
REQUIRED = set(FIXED) | {'script_sha256', 'family_sha256', 'source_point_sha256', 'initial_candidate_sha256', 'checker_sha256', 'runtime_versions', 'dependency_sha256'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def versions():
    return {'python': '.'.join(map(str, sys.version_info[:3])), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'torch': torch.__version__, 'threadpoolctl': threadpoolctl.__version__}


def dependencies():
    return {'python': sha(sys.executable), 'numpy_init': sha(np.__file__),
            'numpy_multiarray': sha(np._core._multiarray_umath.__file__),
            'scipy_init': sha(scipy.__file__), 'torch_init': sha(torch.__file__),
            'torch_C': sha(torch._C.__file__), 'torch_cpu': sha(Path(torch.__file__).parent/'lib/libtorch_cpu.so'),
            'threadpoolctl': sha(threadpoolctl.__file__)}


def guard():
    if not HERE.name.isdigit() or HERE.parent.name != 'runs' or HERE.parent.parent.name != 'experiments':
        raise RuntimeError('execute only from experiments/runs/<reserved-id>')
    root, run_id = HERE.parents[2], int(HERE.name)
    cfg = json.loads((HERE/'config.json').read_text())
    if set(cfg) != REQUIRED:
        raise RuntimeError(f'config-key mismatch: {sorted(set(cfg)^REQUIRED)}')
    for key, expected in FIXED.items():
        if type(cfg[key]) is not type(expected) or cfg[key] != expected:
            raise RuntimeError(f'frozen configuration mismatch: {key}')
    for field, filename in [('script_sha256', 'refine.py'), ('family_sha256', 'family.py'),
                            ('source_point_sha256', cfg['source_point_filename']),
                            ('initial_candidate_sha256', cfg['initial_candidate_filename'])]:
        if cfg[field] != sha(HERE/filename):
            raise RuntimeError(f'frozen input/source hash mismatch: {field}')
    if cfg['checker_sha256'] != sha(root/'check_circuit.py') or cfg['runtime_versions'] != versions() or cfg['dependency_sha256'] != dependencies():
        raise RuntimeError('checker/runtime/dependency mismatch')
    point = json.loads((HERE/cfg['source_point_filename']).read_text())
    if point['parent_run'] != cfg['parent_run'] or point['deletion_ordinal'] != cfg['deletion_ordinal'] or len(point['initial144']) != cfg['parameter_count'] or len(point['topology']) != cfg['expected_cx_count'] or point['candidate_sha256'] != cfg['initial_candidate_sha256']:
        raise RuntimeError('source is not the frozen deletion9 near-hit')
    if point['parent_result_sha256'] != sha(root/'experiments/runs/409/result.json'):
        raise RuntimeError('parent saved result changed')
    with sqlite3.connect(f"file:{root/'experiments.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute('SELECT status,code_sha256,config_json,method,parent_id FROM attempts WHERE id=?', (run_id,)).fetchone()
        parent = db.execute('SELECT status FROM attempts WHERE id=?', (cfg['parent_run'],)).fetchone()
    if row is None or row[0] != 'running' or row[1] != cfg['script_sha256'] or json.loads(row[2]) != cfg or row[3] != cfg['method'] or row[4] != cfg['parent_run']:
        raise RuntimeError('active ledger code/config/method/parent mismatch')
    if parent != ('failed',):
        raise RuntimeError('parent409 is not closedfailed; this must be a fresh run')
    with sqlite3.connect(f"file:{root/'collaboration/board.sqlite3'}?mode=ro", uri=True) as db:
        board = db.execute('SELECT status,owner,run_id FROM hypotheses WHERE id=?', (cfg['board_id'],)).fetchone()
    if board != ('running', 'numerical_sol', run_id) or any(os.environ[key] != '1' for key in THREAD_KEYS):
        raise RuntimeError('board ownership/run link or thread environment mismatch')
    return root, run_id, cfg, point


class Deadline(BaseException):
    pass


def main():
    started = time.monotonic()
    root, run_id, cfg, point = guard()
    with (HERE/'execution_started.json').open('x') as marker:
        json.dump({'run_id': run_id, 'parent_run': cfg['parent_run'], 'config_sha256': sha(HERE/'config.json'), 'started_unix': time.time()}, marker)
    def alarm(signum, frame):
        raise Deadline('internal55-second refinement cap reached')
    signal.signal(signal.SIGALRM, alarm)
    remaining = cfg['internal_seconds']-(time.monotonic()-started)
    if remaining <= 0:
        raise Deadline('preflight exhausted total budget')
    signal.setitimer(signal.ITIMER_REAL, remaining)
    x = np.asarray(point['initial144'], dtype=np.float64)
    best_x = x.copy()
    best_residual = None
    steps = []
    counts = {'residual_calls': 0, 'jacobian_calls': 0, 'svd_calls': 0}
    stop_reason = 'not_started'
    fd_record = None
    try:
        sys.path.insert(0, str(root))
        import check_circuit as checker
        if Path(checker.__file__).resolve() != (root/'check_circuit.py').resolve():
            raise RuntimeError('unexpected shared checker module')
        spec = importlib.util.spec_from_file_location('frozen_eleven_family', HERE/'family.py')
        family = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(family)
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        with threadpool_limits(limits=1):
            word = family.TorchWord(point['topology'])
            mask = ~torch.eye(16, dtype=torch.bool)
            def residual_tensor(v):
                counts['residual_calls'] += 1
                U = word.unitary(v)
                A = U @ word.target @ U.conj().T
                selected = A[mask]
                return torch.cat((selected.real, selected.imag))
            def residual(v):
                out = residual_tensor(torch.tensor(v, dtype=torch.float64)).detach().numpy()
                if out.shape != (cfg['residual_count'],) or not np.isfinite(out).all():
                    raise FloatingPointError('nonfinite or misshaped residual')
                return out
            def jacobian(v):
                counts['jacobian_calls'] += 1
                variable = torch.tensor(v, dtype=torch.float64, requires_grad=True)
                J = torch.autograd.functional.jacobian(residual_tensor, variable, vectorize=True).detach().numpy()
                if J.shape != (cfg['residual_count'], cfg['parameter_count']) or not np.isfinite(J).all():
                    raise FloatingPointError('nonfinite or misshaped residual Jacobian')
                return J
            def metrics(r):
                return {'squared_residual_norm': float(np.dot(r, r)), 'loss': float(np.dot(r, r)/16),
                        'maxoff': float(np.max(np.hypot(r[:240], r[240:])))}
            current = residual(x)
            best_residual = current.copy()
            J = jacobian(x)
            direction = np.cos(np.arange(1, cfg['parameter_count']+1, dtype=np.float64))
            direction /= np.linalg.norm(direction)
            fd = (residual(x+cfg['fd_step']*direction)-residual(x-cfg['fd_step']*direction))/(2*cfg['fd_step'])
            ad = J @ direction
            thresholds = cfg['fd_absolute_tolerance']+cfg['fd_relative_tolerance']*np.abs(fd)
            fd_record = {'direction': direction.tolist(), 'h': cfg['fd_step'], 'ad': ad.tolist(), 'fd': fd.tolist(),
                         'max_error': float(np.max(np.abs(ad-fd))), 'passed': bool(np.all(np.abs(ad-fd) <= thresholds))}
            write(HERE/'directional_jacobian_check.json', fd_record)
            if not fd_record['passed']:
                raise RuntimeError('source residual Jacobian finite-difference sanity check failed')
            write(HERE/'initial_vector.json', {'initial144': x.tolist(), 'metrics': metrics(current)})
            for step_index in range(cfg['max_steps']):
                if metrics(current)['maxoff'] <= cfg['target_maxoff']:
                    stop_reason = 'target_maxoff1e-13_reached'
                    break
                if step_index:
                    J = jacobian(x)
                left, singular, right_t = np.linalg.svd(J, full_matrices=False)
                counts['svd_calls'] += 1
                cutoff = cfg['svd_rcond']*float(singular[0])
                keep = singular > cutoff
                rank = int(np.count_nonzero(keep))
                record = {'step': step_index+1, 'rank': rank, 'rcond': cfg['svd_rcond'], 'cutoff': cutoff,
                          'singular_values': singular.tolist(), 'before': metrics(current), 'line_search': []}
                if not rank:
                    stop_reason = 'zero_retained_svd_rank'
                    steps.append(record)
                    break
                delta = -(right_t[keep].T @ ((left[:, keep].T @ current)/singular[keep]))
                record['delta_norm'] = float(np.linalg.norm(delta))
                accepted = False
                for line_index in range(cfg['max_line_search']):
                    alpha = 2.0**(-line_index)
                    trial_x = x+alpha*delta
                    trial = residual(trial_x)
                    trial_metrics = metrics(trial)
                    improved = trial_metrics['squared_residual_norm'] < record['before']['squared_residual_norm']
                    record['line_search'].append({'alpha': alpha, 'metrics': trial_metrics, 'accepted': bool(improved)})
                    if improved:
                        x, current = trial_x, trial
                        best_x, best_residual = x.copy(), current.copy()
                        record['accepted_alpha'] = alpha
                        accepted = True
                        write(HERE/'best_candidate.json', family.serialize(best_x, point['topology']))
                        write(HERE/'best_vector_checkpoint.json', {'best144': best_x.tolist(), 'metrics': metrics(best_residual)})
                        break
                record['accepted'] = accepted
                steps.append(record)
                write(HERE/'steps.json', steps)
                if not accepted:
                    stop_reason = 'eight_backtracking_trials_no_strict_decrease'
                    break
            else:
                stop_reason = 'max20_newton_steps_reached'
            best_word = family.serialize(best_x, point['topology'])
            write(HERE/'best_candidate.json', best_word)
            initial_checked = checker.evaluate(json.loads((HERE/cfg['initial_candidate_filename']).read_text()), tolerance=cfg['checker_tolerance'])
            best_checked = checker.evaluate(best_word, tolerance=cfg['checker_tolerance'])
            if initial_checked['cnot_count'] != cfg['expected_cx_count'] or best_checked['cnot_count'] != cfg['expected_cx_count']:
                raise RuntimeError('initial or best CX count differs from eleven')
            result = {'status': 'completed', 'run_id': run_id, 'parent_run': cfg['parent_run'], 'deletion_ordinal': cfg['deletion_ordinal'],
                      'initial144': point['initial144'], 'best144': best_x.tolist(), 'initial_checker': initial_checked,
                      'best_checker': best_checked, 'best_residual_metrics': metrics(best_residual), 'stop_reason': stop_reason,
                      'steps': steps, 'counts': counts, 'fd_check_passed': fd_record['passed'],
                      'script_sha256': cfg['script_sha256'], 'family_sha256': cfg['family_sha256'], 'source_point_sha256': cfg['source_point_sha256'],
                      'initial_candidate_sha256': cfg['initial_candidate_sha256'], 'best_candidate_sha256': sha(HERE/'best_candidate.json'),
                      'config_sha256': sha(HERE/'config.json'), 'runtime_versions': versions(), 'dependency_sha256': cfg['dependency_sha256'],
                      'thread_environment': {key: os.environ[key] for key in THREAD_KEYS}, 'torch_threads': torch.get_num_threads(),
                      'torch_interop_threads': torch.get_num_interop_threads(), 'seconds': time.monotonic()-started,
                      'evidence_scope': cfg['evidence_scope']}
            write(HERE/'result.json', result)
            print(json.dumps({'run_id': run_id, 'steps': len(steps), 'stop_reason': stop_reason,
                              'maxoff': best_checked['off_diagonal_error'], 'valid': best_checked['valid_diagonalizer'],
                              'seconds': result['seconds']}, sort_keys=True))
    except BaseException as error:
        signal.setitimer(signal.ITIMER_REAL, 0)
        # Save best vectors, word and evaluated residuals without any new matrix call.
        write(HERE/'terminal_best_vector.json', {'best144': best_x.tolist(), 'best_residual': None if best_residual is None else best_residual.tolist()})
        if 'family' in locals():
            write(HERE/'best_candidate.json', family.serialize(best_x, point['topology']))
        write(HERE/'failure.json', {'run_id': run_id, 'error_type': type(error).__name__, 'error': str(error),
                                  'stop_reason': stop_reason, 'steps': steps, 'counts': counts, 'seconds': time.monotonic()-started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
