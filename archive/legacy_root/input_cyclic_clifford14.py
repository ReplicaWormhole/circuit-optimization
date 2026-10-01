"""Scan rotation-invariant diagonal Clifford input gauges of exact15.

The three generators S^tensor4, product of nearest-neighbor ring CZs,
and CZ(0,2)CZ(1,3) each commute with the four-wire cyclic shift.
Their 4*2*2 combinations test entangling input symmetries excluded from
the previous uniform one-qubit Clifford sweep. This is a bounded compiler
experiment and makes no lower-bound claim.
"""

import json
from pathlib import Path

from qiskit import QuantumCircuit, transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"
BEST = ROOT / "input_cyclic_clifford14_best.json"
RESULT = ROOT / "input_cyclic_clifford14_result.json"


def gauge(s_power, ring, opposite):
    circuit = QuantumCircuit(4)
    for _ in range(s_power):
        for q in range(4):
            circuit.s(q)
    if ring:
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            circuit.cz(a, b)
    if opposite:
        for a, b in ((0, 2), (1, 3)):
            circuit.cz(a, b)
    return circuit


def main():
    incumbent = to_qiskit(json.loads(SOURCE.read_text()))
    rows = []
    best = None
    for s_power in range(4):
        for ring in (0, 1):
            for opposite in (0, 1):
                source = gauge(s_power, ring, opposite).compose(incumbent)
                optimized = transpile(source, basis_gates=["u", "cx"],
                                      optimization_level=3, seed_transpiler=42)
                candidate = from_qiskit(optimized)
                check = evaluate(candidate)
                row = {"s_power": s_power, "ring": ring,
                       "opposite": opposite,
                       "cnot_count": check["cnot_count"],
                       "off_diagonal_error": check["off_diagonal_error"],
                       "valid": check["valid_diagonalizer"]}
                rows.append(row)
                if check["valid_diagonalizer"] and (
                        best is None or row["cnot_count"] < best["cnot_count"]):
                    best = row
                    BEST.write_text(json.dumps(candidate, indent=2) + "\n")
    RESULT.write_text(json.dumps({"trials": len(rows), "best": best,
                                  "rows": rows}, indent=2) + "\n")
    print(json.dumps({"trials": len(rows), "best": best}, indent=2))


if __name__ == "__main__":
    main()
