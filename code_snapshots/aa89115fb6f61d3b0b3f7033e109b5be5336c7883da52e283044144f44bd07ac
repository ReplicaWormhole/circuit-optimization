"""Bounded joint Qiskit synthesis of an exact direct-orbit eigenbasis.

The target matrix is U_block P, where P is run112's two-shear-plus-linear
coordinate permutation and U_block is run119's exact branch Fourier/Bell
transform. This asks the unitary synthesizer to optimize across every
boundary at once; it is a heuristic compiler test, not a CNOT lower bound.
"""

import json
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit
from topology14_orbit_block_rank import full_block


ROOT = Path(__file__).parent
SEEDS = (0, 42)


def main():
    action = json.loads((ROOT / "topology14_orbit_action_result.json").read_text())
    coordinate_map = action["coordinate_map"]
    permutation = np.zeros((16, 16), dtype=complex)
    for x, y in enumerate(coordinate_map):
        permutation[y, x] = 1
    block = np.array(full_block().evalf(), dtype=complex)
    target = block @ permutation
    shift = np.zeros((16, 16), dtype=complex)
    for x in range(16):
        shift[((x & 1) << 3) | (x >> 1), x] = 1
    diagonal = target @ shift @ target.conj().T
    source_offdiag = float(np.max(np.abs(
        diagonal - np.diag(np.diag(diagonal)))))
    assert source_offdiag < 1e-12
    circuit = QuantumCircuit(4)
    circuit.unitary(target, range(4))
    results = []
    for seed in SEEDS:
        compiled = transpile(circuit, basis_gates=["u", "cx"],
                             optimization_level=3,
                             seed_transpiler=seed)
        candidate = from_qiskit(compiled)
        check = evaluate(candidate)
        record = {"seed": seed,
                  "cnot_count": check["cnot_count"],
                  "off_diagonal_error": check["off_diagonal_error"],
                  "valid_diagonalizer": check["valid_diagonalizer"]}
        results.append(record)
        path = ROOT / f"algebraic14_orbit_joint_seed{seed}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
    result = {"source_matrix_off_diagonal_error": source_offdiag,
              "trials": results,
              "best_cnot_count": min(row["cnot_count"] for row in results)}
    (ROOT / "algebraic14_orbit_joint_unitary_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
