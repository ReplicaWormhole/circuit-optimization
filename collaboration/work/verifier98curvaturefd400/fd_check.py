"""Reserved finite-difference validation; run only after run400 is reserved and linked.
Exactly 10 independent NumPy coordinate-matrix evaluations; no autodiff, full
Hessian, fitting, or assignment at perturbed points.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE = ROOT / 'collaboration/work/verifier98fresh/independent_audit.py'
spec = importlib.util.spec_from_file_location('fd_independent_matrix', BASE)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def label_value(u, labels, v):
    return float(np.linalg.norm(u @ v - np.diag(labels) @ u, 'fro') ** 2 / 16)


def off_value(u, v):
    b = u @ v @ u.conj().T
    off = b - np.diag(np.diag(b))
    return float(np.linalg.norm(off, 'fro') ** 2 / 16)


def main():
    manifest_path = HERE / 'fd_manifest.json'
    frozen_path = HERE / 'fd_config.json'
    if (HERE / 'fd_result.json').exists():
        raise RuntimeError('refusing to repeat reserved FD evaluation')
    frozen = json.loads(frozen_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    try:
        here_run_id = int(HERE.name)
    except ValueError as exc:
        raise SystemExit('checker must be placed in its numeric reserved run directory') from exc
    if frozen.get('run_id') != here_run_id:
        raise SystemExit('frozen run ID does not match checker directory')
    if sha(Path(__file__)) != frozen['script_sha256'] or sha(manifest_path) != frozen['manifest_sha256']:
        raise SystemExit('frozen FD script or manifest hash mismatch')
    if frozen.get('external_seconds') != 30 or frozen.get('threads') != 1 or frozen.get('max_matrix_evaluations') != 10:
        raise SystemExit('FD runtime budget differs from freeze')
    if manifest.get('schema') != 'full98-independent-finite-difference-v1' or manifest.get('total_independent_matrix_evaluations') != frozen['max_matrix_evaluations']:
        raise SystemExit('manifest schema/evaluation count mismatch')
    for relative_path, expected_hash in manifest['source_hashes'].items():
        if sha(ROOT / relative_path) != expected_hash:
            raise SystemExit(f'frozen FD source changed: {relative_path}')
    parent_path = ROOT / 'experiments/runs/398/result.json'
    parent = json.loads(parent_path.read_text())
    if sha(parent_path) != manifest['source_result_sha256']:
        raise SystemExit('run398 source result hash mismatch')
    point = np.asarray(manifest['point'], dtype=float)
    if point.shape != (98,) or hashlib.sha256(np.asarray(point, dtype='<f8').tobytes()).hexdigest() != manifest['point_sha256_f64le']:
        raise SystemExit('frozen point malformed or hash mismatch')
    if not np.array_equal(point, np.asarray(parent['point98'], dtype=float)) or parent['seed'] != manifest['seed']:
        raise SystemExit('frozen point/seed differ from source diagnostic')
    records = {record['name']: record for record in parent['records']}
    labels = np.asarray([complex(real, imag) for real, imag in manifest['fixed_labels_real_imag']], dtype=complex)
    source_labels = np.asarray([complex(real, imag) for real, imag in parent['labels_real_imag']], dtype=complex)
    if not np.array_equal(labels, source_labels):
        raise SystemExit('fixed FD labels differ from run398 frozen base assignment')
    for direction_record in manifest['directions']:
        vector = np.asarray(direction_record['direction'], dtype=float)
        source_record = records[direction_record['objective']]
        stored = np.asarray(source_record['eigenvectors_columns'], dtype=float)[:, 0]
        if vector.shape != (98,) or not np.array_equal(vector, stored):
            raise SystemExit('direction is not frozen source eigenvector column 0')
        if hashlib.sha256(np.asarray(vector, dtype='<f8').tobytes()).hexdigest() != direction_record['direction_sha256_f64le']:
            raise SystemExit('direction hash mismatch')
        if abs(float(direction_record['eigenvalue']) - float(source_record['least_eigenvalue'])) > 1e-15:
            raise SystemExit('direction eigenvalue binding mismatch')
    target = base.right_shift()
    calls = []

    def evaluate(point_value, tag):
        if len(calls) >= 10:
            raise RuntimeError('matrix-evaluation cap exceeded')
        unitary = base.coordinate_matrix(np.asarray(point_value, dtype=float))
        unitary_error = float(np.max(np.abs(unitary.conj().T @ unitary - np.eye(16))))
        if unitary_error > 2e-12:
            raise SystemExit(f'nonunitary independent point {tag}: {unitary_error}')
        calls.append({'index': len(calls) + 1, 'tag': tag, 'matrix_unitarity_max_error': unitary_error})
        return unitary

    center_values = {}
    with threadpool_limits(limits=1):
        center_off = evaluate(point, 'center_offdiagonal')
        center_values['offdiagonal'] = off_value(center_off, target)
        center_label = evaluate(point, 'center_fixed_label')
        center_values['fixed_label_branch'] = label_value(center_label, labels, target)
        samples = {}
        for direction_record in manifest['directions']:
            name = direction_record['objective']
            vector = np.asarray(direction_record['direction'], dtype=float)
            for step in manifest['steps']:
                values = {}
                for sign in (-1, 1):
                    tag = f'{name}_h{step:g}_' + ('minus' if sign < 0 else 'plus')
                    unitary = evaluate(point + sign * float(step) * vector, tag)
                    values[sign] = off_value(unitary, target) if name == 'offdiagonal' else label_value(unitary, labels, target)
                center = center_values[name]
                minus, plus = values[-1], values[1]
                first = (plus - minus) / (2 * float(step))
                second = (plus - 2 * center + minus) / (float(step) ** 2)
                gradient = np.asarray(records[name]['gradient'], dtype=float)
                predicted_first = float(np.dot(gradient, vector))
                eigenvalue = float(records[name]['least_eigenvalue'])
                samples[f'{name}|h={step:g}'] = {
                    'direction_source': name,
                    'objective': name,
                    'step': float(step),
                    'minus_value': minus,
                    'center_value': center,
                    'plus_value': plus,
                    'central_first_difference': first,
                    'autodiff_gradient_dot_direction': predicted_first,
                    'first_difference_abs_error': abs(first - predicted_first),
                    'central_second_difference': second,
                    'source_least_eigenvalue': eigenvalue,
                    'second_difference_abs_error_vs_source_eigenvalue': abs(second - eigenvalue),
                }
    if len(calls) != 10:
        raise SystemExit(f'expected exactly 10 independent matrix evaluations, got {len(calls)}')
    output = {
        'schema': manifest['schema'],
        'run_id': here_run_id,
        'frozen_config_sha256': sha(frozen_path),
        'manifest_sha256': sha(manifest_path),
        'fd_script_sha256': sha(Path(__file__)),
        'source_run398_sha256': sha(parent_path),
        'source_point_sha256_f64le': manifest['point_sha256_f64le'],
        'matrix_evaluation_count': len(calls),
        'matrix_evaluations': calls,
        'center_values': center_values,
        'finite_difference_samples': samples,
        'runtime_budget': {
            'external_seconds': frozen['external_seconds'],
            'threads': frozen['threads'],
            'max_matrix_evaluations': frozen['max_matrix_evaluations'],
        },
        'scope': 'Independent finite differences along two frozen least-eigenvector columns for their matching objectives only; no AD/full Hessian/optimizer/assignment solve',
    }
    (HERE / 'fd_result.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
