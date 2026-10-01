"""Test single-wire Schmidt constraints for cross-butterfly output gauges.

The final two lifted Fourier butterflies act on disjoint pairs (0,1) and
(2,3). Left multiplication by any diagonal gate preserves V4
diagonalization. For each cross-pair RZZ gate and delta=k*pi/16, compute
the 4x64 operator-Schmidt matricization across each one-wire cut. Rank 4
on every wire forces at least two incident CNOTs per wire and hence at
least four CNOTs in the whole four-wire block.
"""

import json
import math
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator

from check_circuit import evaluate
from qiskit_crosscheck import to_qiskit


ROOT = Path(__file__).parent
PAIRS = ((0, 2), (0, 3), (1, 2), (1, 3))


def single_wire_rank(unitary, qubit):
    tensor = unitary.reshape((2,) * 8)
    chosen = (qubit, 4 + qubit)
    other = tuple(index for index in range(8) if index not in chosen)
    matrix = tensor.transpose(chosen + other).reshape(4, 64)
    return int(np.linalg.matrix_rank(matrix, tol=1e-9))


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    exact15 = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    assert exact15["gates"][-20:] == baseline["gates"][42:62]
    prefix = exact15["gates"][:-20]
    butterflies = baseline["gates"][42:62]
    results = []
    for pair in PAIRS:
        for k in range(17):
            delta = k * math.pi / 16
            gauge = [] if delta == 0 else [
                {"gate": "cx", "control": pair[0], "target": pair[1]},
                {"gate": "rz", "qubit": pair[1], "theta": delta},
                {"gate": "cx", "control": pair[0], "target": pair[1]},
            ]
            block = {"n": 4, "gates": butterflies + gauge}
            unitary = Operator(to_qiskit(block)).data
            ranks = [single_wire_rank(unitary, q) for q in range(4)]
            check = evaluate({"n": 4, "gates": prefix + block["gates"]})
            results.append({"pair": pair, "k": k,
                            "single_wire_schmidt_ranks": ranks,
                            "full_off_diagonal_error": check["off_diagonal_error"],
                            "full_valid": check["valid_diagonalizer"]})
    print(json.dumps({"cases": len(results),
                      "all_ranks_four": all(r["single_wire_schmidt_ranks"]
                                            == [4, 4, 4, 4] for r in results),
                      "all_full_valid": all(r["full_valid"] for r in results),
                      "results": results}, indent=2))


if __name__ == "__main__":
    main()
