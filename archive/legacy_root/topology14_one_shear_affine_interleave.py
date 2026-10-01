"""Exact controlled-affine gate list and bounded joint Fourier synthesis.

First verify the run136 shared normal form as an exact Qiskit gate circuit
for the one-shear conjugated V4 permutation. Then build the exact 4-qubit
orbit-wise Fourier target, both with and without algebraic run129's RCCX
phase correction. Test direct P+B and interleaved P+L+(B L^-1) placements
under Qiskit level-3 u/cx synthesis. Synthesized angles are numerical;
the target matrices and the action gate list are exact.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator

from algebraic14_one_shear_cost import FORM
from algebraic14_orbit_relative_shears import emit
from algebraic14_orbit_shear_cx import linear_bfs, shortest_form
from check_circuit import right_shift
from topology14_one_shear_orbit_fourier_rank import cycles
from topology14_one_shear_row_gauge_scan import basis_matrix


ROOT = Path(__file__).resolve().parent


def normal_form_action():
    circuit = QuantumCircuit(4)
    circuit.cx(1, 0)  # y3 <- y3 xor y2 = z
    circuit.cx(2, 3)  # L: u1 -> u0
    circuit.cx(1, 2)  # L: u2 -> u1
    circuit.cx(2, 3)  # z0: u1 -> u0
    circuit.cx(3, 1)  # z0: u0 -> u2
    circuit.ccx(0, 3, 1)  # z-controlled correction u0 -> u2
    circuit.cx(0, 3)  # z-controlled X on u0
    circuit.cx(0, 2)  # z-controlled X on u1
    circuit.cx(1, 2)  # L^-1
    circuit.cx(2, 3)
    circuit.cx(1, 0)  # restore y3 = z xor new u2
    return circuit


def coordinate_preprocessor(relative, conjugator):
    circuit = QuantumCircuit(4)
    emit(circuit, FORM, conjugator, relative)
    circuit.cx(3, 1)  # physical 0 -> 2
    circuit.cx(2, 1)  # physical 1 -> 2
    return circuit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    action = json.loads((ROOT / "topology14_one_shear_action_result.json").read_text())
    expected = action["induced_permutation"]
    controlled = normal_form_action()
    action_u = Operator(controlled).data
    action_error = float(np.max(np.abs(action_u - np.eye(16)[:, expected])))
    assert action_error < 1e-12, action_error
    action_compiled = transpile(controlled, basis_gates=["u", "cx"],
                                optimization_level=3, seed_transpiler=args.seed)
    compiled_action_error = float(np.max(np.abs(Operator(action_compiled).data - action_u)))

    order, parents = linear_bfs()
    conjugator = shortest_form(FORM, order, parents)
    orbits = cycles(expected)
    coordinate = action["coordinate_map"]
    defect = json.loads((ROOT / "algebraic14_one_shear_cost_result.json").read_text())[
        "records"][1]["phase_exponents_mod8"]
    phase_corrector = np.array([
        np.exp(-1j * math.pi * defect[coordinate.index(y)] / 4)
        for y in range(16)])
    l = QuantumCircuit(4)
    l.cx(2, 3)
    l.cx(1, 2)
    lmat = Operator(l).data
    records = []
    for relative in (False, True):
        prefix = coordinate_preprocessor(relative, conjugator)
        phase = phase_corrector if relative else np.ones(16, dtype=complex)
        block = basis_matrix(orbits, [list(range(len(orbit))) for orbit in orbits],
                             phase)
        for interleaved in (False, True):
            circuit = prefix.copy()
            target_block = block
            if interleaved:
                circuit.compose(l, inplace=True)
                target_block = block @ lmat.conj().T
            circuit.append(UnitaryGate(target_block), range(4))
            u = Operator(circuit).data
            diagonalized = u @ right_shift(4) @ u.conj().T
            offdiag = float(np.max(np.abs(diagonalized - np.diag(np.diag(diagonalized)))))
            assert offdiag < 1e-9, offdiag
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=args.seed)
            out = Operator(optimized).data
            optimized_diagonalized = out @ right_shift(4) @ out.conj().T
            optimized_offdiag = float(np.max(np.abs(
                optimized_diagonalized - np.diag(np.diag(optimized_diagonalized)))))
            records.append({"relative_phase_shear": relative,
                            "shared_linear_interleaved": interleaved,
                            "prefix_cnot": (11 if relative else 14),
                            "compiled_total_cnot": int(optimized.count_ops().get("cx", 0)),
                            "target_offdiag": offdiag,
                            "compiled_offdiag": optimized_offdiag})
    result = {"controlled_affine_action_gate_list": [
                  {"name": instruction.operation.name,
                   "qubits_qiskit_little_endian": [controlled.find_bit(q).index
                                                   for q in instruction.qubits]}
                  for instruction in controlled.data],
              "action_permutation_error": action_error,
              "compiled_action_cnot": int(action_compiled.count_ops().get("cx", 0)),
              "compiled_action_error": compiled_action_error,
              "joint_fourier_records": records,
              "synthesis_basis": ["u", "cx"],
              "exact_status": "controlled affine gate list exact; Fourier target exact cyclotomic matrix; transpiled u angles numerical"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"action_error": action_error,
                      "compiled_action_cnot": result["compiled_action_cnot"],
                      "compiled_action_error": compiled_action_error,
                      "joint_fourier_records": records}, indent=2))


if __name__ == "__main__":
    main()
