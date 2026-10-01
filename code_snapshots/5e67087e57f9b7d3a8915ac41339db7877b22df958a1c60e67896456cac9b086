"""One reserved full98 curvature diagnostic, no optimizer."""
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
from scipy.optimize import linear_sum_assignment
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    out = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out


def main():
    cfg = json.loads((HERE / 'config.json').read_text())
    if (HERE / 'result.json').exists():
        raise RuntimeError('Refusing to repeat a reserved diagnostic')
    assert digest(__file__) == cfg['code_sha256']
    for p, h in cfg['dependencies'].items():
        assert digest(ROOT / p) == h, p
    assert scipy.__version__ == cfg['scipy_version']
    family = module('curvature_family', ROOT / cfg['family_module'])
    assignment_module = module('curvature_assignment', ROOT / cfg['assignment_module'])
    parent = json.loads((ROOT / cfg['parent_result']).read_text())
    row = next(r for r in parent['rows'] if r['seed'] == cfg['seed'])
    point = np.asarray(row['best'], dtype=float)
    assert point.shape == (98,) and np.isfinite(point).all()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    started = time.monotonic()
    deadline = started + cfg['internal_seconds']
    records = []
    termination = 'all_diagnostics_finished'
    class BudgetStop(Exception):
        pass
    def budget():
        if time.monotonic() >= deadline:
            raise BudgetStop()

    with threadpool_limits(limits=1):
        u = family.matrix(torch.tensor(point, dtype=torch.float64)).detach().numpy()
        v_np = family.right_shift(4).astype(complex)
        v = torch.tensor(v_np, dtype=torch.complex128)
        labels, columns, label_cost = assignment_module.assignment(u, v_np)
        uv = u @ v_np
        slots = assignment_module.SLOTS
        costs = np.sum(np.abs(uv[:, None, :] - slots[None, :, None] * u[:, None, :]) ** 2, axis=2)
        alternatives = []
        gap = None
        d = torch.tensor(labels, dtype=torch.complex128)
        def off(z):
            budget()
            return family.objective(z)
        def branch(z):
            budget()
            m = family.matrix(z)
            residual = m @ v - d[:, None] * m
            return (residual.abs() ** 2).sum().real / 16
        try:
            for i in range(16):
                budget()
                constrained = costs.copy()
                constrained[i, slots == labels[i]] = np.inf
                rows, cols = linear_sum_assignment(constrained)
                alternatives.append({'banned_row': i, 'columns': cols.tolist(),
                                     'normalized_cost': float(constrained[rows, cols].sum() / 16)})
            gap = min(r['normalized_cost'] for r in alternatives) - label_cost
            for name, fun in [('offdiagonal', off), ('fixed_label_branch', branch)]:
                budget()
                t0 = time.monotonic()
                x = torch.tensor(point, dtype=torch.float64, requires_grad=True)
                value = fun(x)
                gradient_tensor = torch.autograd.grad(value, x, create_graph=True)[0]
                gradient = gradient_tensor.detach().numpy()
                rec = {'name': name, 'loss': float(value.detach()), 'gradient': gradient.tolist(),
                       'gradient_norm': float(np.linalg.norm(gradient)), 'complete': False,
                       'partial_hessian_rows': [], 'hessian_rows_completed': 0}
                records.append(rec)
                for j in range(98):
                    budget()
                    hrow = torch.autograd.grad(gradient_tensor[j], x, retain_graph=True)[0].detach().numpy()
                    assert np.isfinite(hrow).all()
                    rec['partial_hessian_rows'].append(hrow.tolist())
                    rec['hessian_rows_completed'] += 1
                budget()
                h = np.asarray(rec['partial_hessian_rows'])
                hs = (h + h.T) / 2
                eig, vec = np.linalg.eigh(hs)
                for j in range(98):
                    i = int(np.argmax(np.abs(vec[:, j])))
                    if vec[i, j] < 0:
                        vec[:, j] *= -1
                rec.update(hessian=h.tolist(), hessian_symmetry_max=float(np.max(np.abs(h-h.T))),
                           eigenvalues=eig.tolist(), eigenvectors_columns=vec.tolist(),
                           least_eigenvalue=float(eig[0]),
                           least_eigenpair_residual=float(np.linalg.norm(hs @ vec[:, 0] - eig[0] * vec[:, 0])),
                           negative_trigger=bool(eig[0] < cfg['negative_trigger']),
                           elapsed_seconds=time.monotonic()-t0, complete=True)
                rec.pop('partial_hessian_rows')
        except BudgetStop:
            termination = 'internal_time_limit'
        native = family.serialize(point)
        compiled = family.helper.compile_native(native)
        for name, payload in [('base_native', native), ('base_compiled', compiled)]:
            (HERE / (name + '.json')).write_text(json.dumps(payload, indent=2) + '\n')
        out = {'run_id': int(HERE.name), 'config': cfg, 'point98': point.tolist(),
               'point_sha256_f64le': digest_bytes(point), 'seed': cfg['seed'],
               'labels_real_imag': [[float(z.real), float(z.imag)] for z in labels],
               'columns': columns.tolist(), 'label_cost': label_cost,
               'distinct_root_assignment_gap': gap, 'assignment_alternatives': alternatives,
               'records': records, 'termination': termination,
               'elapsed_seconds': time.monotonic()-started,
               'base_native_sha256': digest(HERE/'base_native.json'),
               'base_compiled_sha256': digest(HERE/'base_compiled.json'),
               'base_checker': family.evaluate(compiled),
               'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'torch': torch.__version__},
               'scope': 'Floating-point curvature at one fixed full98 point; no fit, local-minimum proof, exclusion or certificate'}
        (HERE / 'result.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'run_id': out['run_id'], 'gap': gap, 'elapsed_seconds': out['elapsed_seconds'],
                      'termination': termination,
                      'diagnostics': [{k: r[k] for k in ('name', 'loss', 'gradient_norm', 'complete', 'least_eigenvalue', 'negative_trigger') if k in r} for r in records]}, indent=2))


def digest_bytes(point):
    return hashlib.sha256(np.asarray(point, dtype='<f8').tobytes()).hexdigest()


if __name__ == '__main__':
    main()
