"""Finite objective fixtures for the dynamic-eigenlabel full98 run.

No optimizer or derivatives. Uses the preserved independent NumPy/SciPy
coordinate/gate evaluator and raw pairwise row costs for assignment.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[3]
BASE_PATH = ROOT / "collaboration/work/verifier98fresh/independent_audit.py"
BASE_SPEC = importlib.util.spec_from_file_location("verifier98_base_for_label", BASE_PATH)
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)

LABELS = np.repeat(np.array([1, -1, 1j, -1j], dtype=complex), [6, 4, 3, 3])


def raw_cost_matrix(unitary: np.ndarray, target: np.ndarray | None = None) -> np.ndarray:
    """Frozen-script raw costs ||(UV)_i-lambda_j U_i||^2, one row per i."""
    v = base.right_shift() if target is None else target
    uv = unitary @ v
    return np.sum(np.abs(uv[:, None, :] - LABELS[None, :, None] * unitary[:, None, :]) ** 2, axis=2)


def assign_labels(unitary: np.ndarray, target: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    costs = raw_cost_matrix(unitary, target)
    row_ind, col_ind = linear_sum_assignment(costs)
    d = np.empty(16, dtype=complex)
    d[row_ind] = LABELS[col_ind]
    return d, costs, np.stack((row_ind, col_ind))


def objective_components(unitary: np.ndarray) -> dict:
    v = base.right_shift()
    conjugated = unitary @ v @ unitary.conj().T
    d, costs, indices = assign_labels(unitary, v)
    D = np.diag(d)
    direct = unitary @ v - D @ unitary
    conjugated_residual = conjugated - D
    off = conjugated - np.diag(np.diag(conjugated))
    label_direct = float(np.linalg.norm(direct, "fro") ** 2 / 16)
    label_conjugated = float(np.linalg.norm(conjugated_residual, "fro") ** 2 / 16)
    raw_assigned_cost = float(costs[indices[0], indices[1]].sum() / 16)
    return {
        "labels_by_row_real_imag": [[float(z.real), float(z.imag)] for z in d],
        "assignment_rows": indices[0].tolist(),
        "assignment_columns": indices[1].tolist(),
        "label_residual_direct": label_direct,
        "label_residual_conjugated": label_conjugated,
        "raw_assignment_cost_normalized": raw_assigned_cost,
        "original_offdiagonal_loss": float(np.linalg.norm(off, "fro") ** 2 / 16),
        "original_offdiagonal_max": float(np.max(np.abs(off))),
        "identity_abs_difference": abs(label_direct - label_conjugated),
        "cost_abs_difference": abs(label_direct - raw_assigned_cost),
        "raw_cost_sha256_f64le": hashlib.sha256(np.asarray(costs, dtype="<f8").tobytes()).hexdigest(),
    }


def haar_fixture(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16))
    q, r = np.linalg.qr(z)
    phases = np.diag(r)
    q = q @ np.diag(np.conjugate(phases / np.abs(phases)))
    return q


def multiplicity_audit() -> dict:
    values = np.linalg.eigvals(base.right_shift())
    targets = np.array([1, -1, 1j, -1j], dtype=complex)
    counts = [int(np.sum(np.abs(values - target) < 1e-10)) for target in targets]
    if counts != [6, 4, 3, 3]:
        raise AssertionError(f"unexpected cycle eigenvalue multiplicities: {counts}")
    return {"eigenvalue_targets_real_imag": [[float(z.real), float(z.imag)] for z in targets], "multiplicities": counts,
            "sum": sum(counts), "eigenvalue_match_max_error": float(max(min(abs(v-target) for target in targets) for v in values))}


def fixtures() -> dict:
    mult = multiplicity_audit()
    # Identity creates assignment ties; exact repeatability is checked on the
    # same raw cost array and pinned SciPy implementation.
    identity = np.eye(16, dtype=complex)
    tie_cost = raw_cost_matrix(identity, base.right_shift())
    tie1 = linear_sum_assignment(tie_cost)
    tie2 = linear_sum_assignment(tie_cost.copy())
    if not (np.array_equal(tie1[0], tie2[0]) and np.array_equal(tie1[1], tie2[1])):
        raise AssertionError("SciPy tie fixture assignment is not repeatable")
    tie_digest = hashlib.sha256(np.asarray(tie_cost, dtype="<f8").tobytes()).hexdigest()

    random_u = haar_fixture(926162099)
    random_report = objective_components(random_u)
    eigvals, eigvecs = np.linalg.eig(base.right_shift())
    # V is normal; its eigenvectors form an orthonormal basis up to roundoff.
    eig_u = eigvecs.conj().T
    exact_report = objective_components(eig_u)
    unitarity = float(np.max(np.abs(eig_u.conj().T @ eig_u - np.eye(16))))
    if unitarity > 1e-12:
        raise AssertionError(f"eigenbasis fixture is not unitary: {unitarity}")
    if exact_report["label_residual_direct"] > 1e-12:
        raise AssertionError("eigenbasis fixture failed to reach a zero labeled residual")
    if random_report["identity_abs_difference"] > 1e-12 or random_report["cost_abs_difference"] > 1e-12:
        raise AssertionError("random-unitary identity/cost fixture mismatch")
    return {
        "scipy_version": scipy.__version__,
        "multiplicity_audit": mult,
        "tie_fixture": {"cost_sha256_f64le": tie_digest,
                        "assignment_rows": tie1[0].tolist(),
                        "assignment_columns": tie1[1].tolist(),
                        "repeatable_on_identical_input": True},
        "random_unitary_fixture": random_report,
        "eigenbasis_fixture": {**exact_report, "unitarity_max_error": unitarity},
        "scope": "finite objective-identity, multiplicity, Hungarian raw-cost, and deterministic tie fixtures; no fit or derivatives",
    }


if __name__ == "__main__":
    print(json.dumps(fixtures(), indent=2))
