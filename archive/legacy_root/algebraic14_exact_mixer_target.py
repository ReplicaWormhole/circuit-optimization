"""Exact simple output-gauge target behind the numerical 14-CX circuit.

Let U15 be the verified exact15 diagonalizer. Define M by two disjoint
2x2 blocks on source/output rows with equal V4 eigenvalues:
  source columns (6,9) -> output rows (5,10),
  source columns (5,10) -> output rows (6,9),
each block [[sqrt(3)/2,-i/2],[1/2,i*sqrt(3)/2]].
M is exactly unitary and maps D15 to a diagonal D14, so M U15 exactly
diagonalizes V4. Compare the numerical 14-CX candidate to M U15 up to
free output-row phases; this does not exactify its gate angles.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp
from qiskit.quantum_info import Operator

from exact_check import exact_eigenvalue_labels
from qiskit_crosscheck import to_qiskit


ROOT = Path(__file__).parent


def exact_mixer():
    c, s = sp.sqrt(3) / 2, sp.Rational(1, 2)
    block = sp.Matrix([[c, -sp.I * s], [s, sp.I * c]])
    m = sp.eye(16)
    for rows, cols in (((5, 10), (6, 9)), ((6, 9), (5, 10))):
        for r in rows:
            m[r, r] = 0
        for a, row in enumerate(rows):
            for b, column in enumerate(cols):
                m[row, column] = block[a, b]
    return m


def main():
    exact15 = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    numerical14 = json.loads(
        (ROOT / "algebraic14_mixer_p75_gn_candidate.json").read_text())
    certificate = exact_eigenvalue_labels(exact15)
    assert certificate["exact_diagonalizer"]
    labels15 = certificate["output_labels"]
    roots = (1, sp.I, -1, -sp.I)
    d15 = sp.diag(*(roots[label] for label in labels15))
    m = exact_mixer()
    unitary_exact = all(
        sp.simplify(value) == 0 for value in m * m.conjugate().T - sp.eye(16))
    assert unitary_exact
    d14 = m * d15 * m.conjugate().T
    diagonal_exact = all(
        sp.simplify(d14[i, j]) == 0
        for i in range(16) for j in range(16) if i != j)
    assert diagonal_exact
    labels14 = [
        next(label for label, root in enumerate(roots)
             if sp.simplify(d14[j, j] - root) == 0)
        for j in range(16)]
    u15 = Operator(to_qiskit(exact15)).data
    u14 = Operator(to_qiskit(numerical14)).data
    m_np = np.array(m.evalf(), dtype=complex)
    exact_target_numeric = m_np @ u15
    overlap = u14 @ exact_target_numeric.conj().T
    row_phases = np.diag(overlap) / np.abs(np.diag(overlap))
    adjusted_error = float(np.max(np.abs(
        u14 - row_phases[:, None] * exact_target_numeric)))
    overlap_offdiag = float(np.max(np.abs(
        overlap - np.diag(np.diag(overlap)))))
    result = {
        "exact_mixer_unitary": unitary_exact,
        "exact_conjugated_eigenvalue_matrix_diagonal": diagonal_exact,
        "source_labels": labels15,
        "target_labels": labels14,
        "mixer_block": [["sqrt(3)/2", "-I/2"],
                        ["1/2", "I*sqrt(3)/2"]],
        "numerical14_vs_exact_target_row_phase_error": adjusted_error,
        "output_overlap_offdiag_error": overlap_offdiag,
        "scope": "Exact target matrix, numerical 14-CX implementation",
    }
    (ROOT / "algebraic14_exact_mixer_target_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
