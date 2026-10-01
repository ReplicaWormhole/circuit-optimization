"""Hand-structured exact parity-sector Fourier block after shared affine L.

In L coordinates, odd-parity (z=1) action is decrement on (u0,u1) with
spectator u2. Even-parity (z=0) action has decrement on (u2,u0) when u1=1
and CX(u0,u2) when u1=0. The exact block applies one controlled two-bit
QFT in each four-cycle branch and one controlled H for the two-cycle branch.
Qiskit compiles those exact high-level controlled gates to u/cx numerically.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator

from algebraic14_one_shear_cost import FORM
from algebraic14_orbit_relative_shears import emit
from algebraic14_orbit_shear_cx import linear_bfs, shortest_form
from check_circuit import right_shift
from topology14_orbit_branch_basis import CP_PLUS, H0, H1


ROOT = Path(__file__).resolve().parent


def sector_fourier():
    qft = UnitaryGate(np.array(H1 * CP_PLUS * H0, dtype=complex), label="QFT2")
    h = UnitaryGate(np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2),
                    label="H")
    circuit = QuantumCircuit(4)
    circuit.cx(1, 0)  # extract z = y2 xor y3 onto physical q0
    circuit.cx(2, 3)  # shared L: u1 -> u0
    circuit.cx(1, 2)  # shared L: u2 -> u1

    # z=1: QFT on ordered (u0,u1) = physical (q3,q2).
    circuit.append(qft.control(1), [0, 2, 3])

    # z=0: open parity control, then branch on u1=physical q2.
    circuit.x(0)
    circuit.x(2)
    circuit.append(h.control(2), [0, 2, 1])  # u1=0: H on u2
    circuit.x(2)
    circuit.append(qft.control(2), [0, 2, 3, 1])  # u1=1: QFT on (u2,u0)
    circuit.x(0)
    return circuit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    order, parents = linear_bfs()
    conjugator = shortest_form(FORM, order, parents)
    prefix = QuantumCircuit(4)
    emit(prefix, FORM, conjugator, False)
    prefix.cx(3, 1)
    prefix.cx(2, 1)
    block = sector_fourier()
    full = prefix.compose(block)
    u = Operator(full).data
    diagonalized = u @ right_shift(4) @ u.conj().T
    offdiag = float(np.max(np.abs(diagonalized - np.diag(np.diag(diagonalized)))))
    assert offdiag < 1e-9, offdiag
    compiled_block = transpile(block, basis_gates=["u", "cx"],
                               optimization_level=3, seed_transpiler=args.seed)
    compiled_full = transpile(full, basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=args.seed)
    compiled_u = Operator(compiled_full).data
    compiled_diag = compiled_u @ right_shift(4) @ compiled_u.conj().T
    compiled_offdiag = float(np.max(np.abs(
        compiled_diag - np.diag(np.diag(compiled_diag)))))
    result = {"exact_high_level_gate_sequence": [
                  {"name": item.operation.name,
                   "qubits_qiskit_little_endian": [full.find_bit(q).index
                                                   for q in item.qubits]}
                  for item in full.data],
              "prefix_cnot": 14,
              "structured_block_compiled_cnot": int(compiled_block.count_ops().get("cx", 0)),
              "structured_full_compiled_cnot": int(compiled_full.count_ops().get("cx", 0)),
              "offdiag": offdiag, "compiled_offdiag": compiled_offdiag,
              "exact_status": "abstract controlled QFT/H and CCX gates exact; u/cx compilation numerical"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "exact_high_level_gate_sequence"}, indent=2))


if __name__ == "__main__":
    main()
