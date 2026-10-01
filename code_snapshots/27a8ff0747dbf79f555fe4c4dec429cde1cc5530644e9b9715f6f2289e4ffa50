"""One reserved full98 eigenvalue-assignment fit; numerical evidence only."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path
import numpy as np
import scipy
import torch
from scipy.optimize import linear_sum_assignment, minimize
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
SLOTS = np.array([1] * 6 + [-1] * 4 + [1j] * 3 + [-1j] * 3, dtype=complex)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assignment(u, v):
    """Unperturbed exact-spectrum row costs; pinned SciPy resolves tied optima."""
    uv = u @ v
    cost = np.sum(np.abs(uv[:, None, :] - SLOTS[None, :, None] * u[:, None, :]) ** 2, axis=2)
    rows, columns = linear_sum_assignment(cost)
    assert np.array_equal(rows, np.arange(16))
    labels = SLOTS[columns]
    return labels, columns, float(cost[rows, columns].sum() / 16)


def main():
    cfg = json.loads((HERE / 'config.json').read_text())
    if (HERE / 'result.json').exists():
        raise RuntimeError('Refusing to repeat a reserved fit')
    assert digest(__file__) == cfg['code_sha256']
    for path, sha in cfg['dependencies'].items():
        assert digest(ROOT / path) == sha, path
    assert scipy.__version__ == cfg['scipy_version']
    parent = json.loads((ROOT / cfg['parent_result']).read_text())
    row = next(r for r in parent['rows'] if r['seed'] == cfg['seed'])
    initial = np.asarray(row['best'], dtype=float)
    assert initial.shape == (98,) and np.isfinite(initial).all()
    spec = importlib.util.spec_from_file_location('full98label', ROOT / cfg['family_module'])
    family = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(family)
    v_np = family.right_shift(4).astype(complex)
    v = torch.tensor(v_np, dtype=torch.complex128)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    start = time.monotonic()
    deadline = start + cfg['internal_seconds']
    calls = 0
    history = []
    best = None
    best_original = None
    initial_record = None
    last = None
    opt = None
    termination = 'optimizer_returned'

    class BudgetStop(Exception):
        pass

    class NearHit(Exception):
        pass

    def fun(z):
        nonlocal calls, best, best_original, last, initial_record
        if time.monotonic() >= deadline or calls >= cfg['maxfun']:
            raise BudgetStop()
        x = torch.tensor(z, dtype=torch.float64, requires_grad=True)
        u = family.matrix(x)
        u_np = u.detach().numpy()
        labels, columns, assigned_cost = assignment(u_np, v_np)
        d = torch.tensor(labels, dtype=torch.complex128)
        residual = u @ v - d[:, None] * u
        value = (residual.abs() ** 2).sum().real / 16
        gradient = torch.autograd.grad(value, x)[0].detach().numpy()
        loss = float(value.detach())
        b = u_np @ v_np @ u_np.conj().T
        original = float(np.sum(np.abs(b - np.diag(np.diag(b))) ** 2) / 16)
        if not np.isfinite(loss) or not np.isfinite(gradient).all():
            raise FloatingPointError('Nonfinite evaluation')
        assert abs(loss - assigned_cost) < 1e-12
        calls += 1
        last = {'point': np.asarray(z).tolist(), 'label_loss': loss, 'original_loss': original,
                'columns': columns.tolist(), 'gradient_norm': float(np.linalg.norm(gradient)),
                'evaluation': calls}
        history.append({k: last[k] for k in ('label_loss', 'original_loss', 'columns', 'gradient_norm', 'evaluation')})
        if initial_record is None:
            initial_record = last.copy()
        if best is None or loss < best['label_loss']:
            best = last.copy()
        if best_original is None or original < best_original['original_loss']:
            best_original = last.copy()
        if loss < cfg['nearhit_loss']:
            raise NearHit()
        return loss, gradient

    with threadpool_limits(limits=1):
        try:
            opt = minimize(fun, initial, jac=True, method='L-BFGS-B', options={
                'maxiter': cfg['maxiter'], 'maxfun': cfg['maxfun'], 'maxls': cfg['maxls'],
                'ftol': cfg['ftol'], 'gtol': cfg['gtol']})
        except BudgetStop:
            termination = 'time_or_evaluation_budget'
        except NearHit:
            termination = 'promising_nearhit_requires_independent_verification'
        except FloatingPointError:
            termination = 'nonfinite_evaluation'
        if best is None:
            raise RuntimeError('No evaluated point; close ledger without a candidate')
        endpoints = {'initial': initial_record, 'best': best, 'best_original': best_original, 'last_evaluated': last}
        if opt is not None:
            endpoints['optimizer_return'] = {'point': opt.x.tolist()}
        for name, record in endpoints.items():
            point = np.asarray(record['point'])
            native = family.serialize(point)
            compiled = family.helper.compile_native(native)
            for kind, payload in [('native', native), ('compiled', compiled)]:
                path = HERE / f'{name}_{kind}.json'
                path.write_text(json.dumps(payload, indent=2) + '\n')
                record[kind] = path.name
                record[kind + '_sha256'] = digest(path)
            record['checker'] = family.evaluate(compiled)
    out = {'run_id': int(HERE.name), 'config': cfg, 'seed': cfg['seed'], 'initial': initial.tolist(),
           'endpoints': endpoints, 'history': history, 'evaluations': calls,
           'nit': None if opt is None else int(opt.nit), 'termination': termination,
           'optimizer_message': None if opt is None else str(opt.message),
           'elapsed_seconds': time.monotonic() - start,
           'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'torch': torch.__version__},
           'scope': 'One bounded piecewise-smooth objective fit; no exact certificate or exclusion'}
    (HERE / 'result.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ('run_id', 'seed', 'evaluations', 'nit', 'termination', 'elapsed_seconds')}, indent=2))
    print(json.dumps({'best_label_loss': best['label_loss'], 'best_original_loss': best['original_loss'],
                      'best_checker': best['checker']}, indent=2))


if __name__ == '__main__':
    main()
