"""Finite-difference Gauss-Newton correction of p=3/4 14-CX mixer.

Augment UV4=DU residual by exact output-gauge constraints
|(U U15†)[5,6]|²=|(U U15†)[6,5]|²=3/4. Start from run155's
checker-valid candidate. This is a numerical certificate of a simple
mixer gauge, not exact gate-angle certification.
"""

import json
from pathlib import Path

import numpy as np

from check_circuit import (cnot, embedded_one_qubit, evaluate,
                           one_qubit_matrix)
from delete14_nearmiss_refine import (
    ROOTS, V, emit, parse, unitary)


ROOT = Path(__file__).parent
TARGET = 0.75
DIFFERENCE_STEP = 1e-5
RCOND = 1e-12
MAXITER = 8


def gate_matrix(gate):
    if gate["gate"] == "cx":
        return cnot(gate["control"], gate["target"], 4)
    return embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4)


def circuit_matrix(candidate):
    u = np.eye(16, dtype=complex)
    for gate in candidate["gates"]:
        u = gate_matrix(gate) @ u
    return u


def main():
    source = json.loads((ROOT / "algebraic14_mixer_p75.json").read_text())
    old = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    u15 = circuit_matrix(old)
    prefix, prefix_u, edges, angles = parse(source)
    cx_matrices = [cnot(*edge, 4) for edge in edges]
    current = unitary(prefix_u, cx_matrices, angles)
    diagonal = np.diag(current @ V @ current.conj().T)
    indices = np.argmin(abs(diagonal[:, None] - ROOTS[None, :]), axis=1)
    labels = ROOTS[indices]

    def residual(flat):
        u = unitary(prefix_u, cx_matrices, flat.reshape(9, 4, 3))
        equation = u @ V - labels[:, None] * u
        m = u @ u15.conj().T
        gauge = np.array([abs(m[5, 6]) ** 2 - TARGET,
                          abs(m[6, 5]) ** 2 - TARGET])
        return np.concatenate((equation.real.ravel(),
                               equation.imag.ravel(), gauge))

    history = []
    for iteration in range(MAXITER):
        flat = angles.ravel()
        r = residual(flat)
        norm = float(np.linalg.norm(r))
        maximum = float(np.max(abs(r)))
        if maximum < 1e-13:
            history.append({"iteration": iteration,
                            "norm": norm, "maximum": maximum,
                            "stopped": "threshold"})
            break
        jac = np.empty((len(r), len(flat)), dtype=float)
        for coordinate in range(len(flat)):
            plus, minus = flat.copy(), flat.copy()
            plus[coordinate] += DIFFERENCE_STEP
            minus[coordinate] -= DIFFERENCE_STEP
            jac[:, coordinate] = (
                residual(plus) - residual(minus)) / (2 * DIFFERENCE_STEP)
        delta, _, rank, singular = np.linalg.lstsq(jac, -r, rcond=RCOND)
        accepted = False
        for factor in (1.0, 0.5, 0.25, 0.125):
            trial = flat + factor * delta
            trial_norm = float(np.linalg.norm(residual(trial)))
            if trial_norm < norm:
                angles = trial.reshape(9, 4, 3)
                accepted = True
                break
        history.append({"iteration": iteration, "norm": norm,
                        "maximum": maximum, "rank": int(rank),
                        "step_norm": float(np.linalg.norm(delta)),
                        "accepted_factor": factor if accepted else None,
                        "next_norm": trial_norm if accepted else None})
        if not accepted:
            break
    candidate = emit(prefix, edges, angles)
    path = ROOT / "algebraic14_mixer_p75_gn_candidate.json"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    u = circuit_matrix(candidate)
    m = u @ u15.conj().T
    result = {"target_probability": TARGET,
              "achieved_probabilities": [float(abs(m[5, 6]) ** 2),
                                         float(abs(m[6, 5]) ** 2)],
              "fixed_label_max_residual": float(np.max(np.abs(
                  u @ V - labels[:, None] * u))),
              "offdiag_error": check["off_diagonal_error"],
              "valid": check["valid_diagonalizer"],
              "cnot_count": check["cnot_count"],
              "history": history,
              "candidate": str(path)}
    (ROOT / "algebraic14_mixer_p75_gn_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
