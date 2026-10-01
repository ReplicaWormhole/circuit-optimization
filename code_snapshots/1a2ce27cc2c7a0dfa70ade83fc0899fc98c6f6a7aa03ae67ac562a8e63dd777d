"""Search single output controlled-phase gauges beyond Clifford CZ.

An output-diagonal CP commutes with the diagonal eigenvalue matrix D, so
it preserves V4 diagonalization. The scan is restricted to one pair and
angles ±pi/8, ±pi/4, ±pi/2; it is not an optimality argument.
"""

import itertools
import json
import math
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
PAIRS = tuple(itertools.combinations(range(4), 2))
ANGLES = (-math.pi/2, -math.pi/4, -math.pi/8,
          math.pi/8, math.pi/4, math.pi/2)


def main():
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    rows = []
    best = None
    for pair in PAIRS:
        for angle in ANGLES:
            circuit = to_qiskit(base)
            circuit.cp(angle, 3 - pair[0], 3 - pair[1])
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = from_qiskit(optimized)
            check = evaluate(candidate)
            row = {"pair": pair, "angle_over_pi": angle / math.pi,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid": check["valid_diagonalizer"]}
            rows.append(row)
            if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                    best["cnot_count"], best["off_diagonal_error"]):
                best = row
                (ROOT / "output_cp_gauge14_best.json").write_text(
                    json.dumps(candidate, indent=2) + "\n")
    result = {"trials": len(rows),
              "minimum_cnot": min(row["cnot_count"] for row in rows),
              "best": best, "rows": rows}
    (ROOT / "output_cp_gauge14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "minimum_cnot", "best")}, indent=2))


if __name__ == "__main__":
    main()
