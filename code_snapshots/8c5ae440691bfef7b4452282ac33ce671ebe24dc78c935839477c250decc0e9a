"""Compare parity-conditioned Fock Fourier factors to the 18-CX baseline.

Determine whether the baseline eight-CX diagonal prefix is the same
parity twist G up to a shift-invariant input gauge, and whether the
ten-CX Fourier tail is a monomial output relabeling of periodic Fock FFT.
This informs which butterfly may absorb an open twist network.
"""

import json
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator

from algebraic14_fermion_twist_factor import periodic_transform, phase_integer
from algebraic14_parity_fermion_fourier import FIELD
from qiskit_crosscheck import to_qiskit


ROOT = Path(__file__).parent
ZETA = np.exp(1j * np.pi / 4)


def numeric_field(matrix):
    return np.array([[complex(FIELD.to_sympy(z).evalf())
                      for z in row] for row in matrix], dtype=complex)


def rotate(x):
    return ((x & 1) << 3) | (x >> 1)


def main():
    gates = json.loads((ROOT / "baseline_18.json").read_text())["gates"]
    prefix = Operator(to_qiskit({"n": 4, "gates": gates[:15]})).data
    tail = Operator(to_qiskit({"n": 4, "gates": gates[15:]})).data
    u0 = numeric_field(periodic_transform())
    g = np.diag([ZETA ** phase_integer(x) for x in range(16)])
    prefix_offdiag = float(np.max(np.abs(
        prefix - np.diag(np.diag(prefix)))))
    relative = prefix @ g.conj().T
    relative_values = np.diag(relative)
    commute_error = max(abs(relative_values[rotate(x)] -
                            relative_values[x]) for x in range(16))
    change = tail @ u0.conj().T
    monomial_rows = [int(np.count_nonzero(abs(row) > 1e-8))
                     for row in change]
    result = {
        "prefix_cx": sum(gate["gate"] == "cx" for gate in gates[:15]),
        "tail_cx": sum(gate["gate"] == "cx" for gate in gates[15:]),
        "prefix_offdiag": prefix_offdiag,
        "prefix_vs_twist_shift_invariant_error": float(commute_error),
        "tail_vs_periodic_fock_overlap_row_nonzero_counts": monomial_rows,
        "tail_vs_periodic_fock_monomial": all(count == 1
                                             for count in monomial_rows)}
    (ROOT / "algebraic14_fermion_baseline_relation_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
