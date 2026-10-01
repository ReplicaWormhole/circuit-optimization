"""Search input dihedral qubit symmetries of the 15-CNOT V4 circuit.

Dihedral wire permutations conjugate the right shift to itself or its inverse.
Precomposing a V4 diagonalizer with one therefore preserves diagonalization,
possibly changing output eigenvalue labels. This is an eight-case scan.
"""

import json
from pathlib import Path

from qiskit import QuantumCircuit, transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent


def swaps_for(target):
    current = list(range(4))
    swaps = []
    for index in range(4):
        if current[index] == target[index]:
            continue
        partner = current.index(target[index])
        current[index], current[partner] = current[partner], current[index]
        swaps.append((index, partner))
    return swaps


def main():
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    rows = []
    best = None
    for sign in (1, -1):
        for offset in range(4):
            permutation = tuple((sign * index + offset) % 4 for index in range(4))
            prefix = QuantumCircuit(4)
            for a, b in swaps_for(permutation):
                prefix.swap(3 - a, 3 - b)
            prefix.compose(to_qiskit(base), inplace=True)
            optimized = transpile(prefix, basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = from_qiskit(optimized)
            check = evaluate(candidate)
            row = {"sign": sign, "offset": offset,
                   "permutation": permutation,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid": check["valid_diagonalizer"]}
            rows.append(row)
            if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                    best["cnot_count"], best["off_diagonal_error"]):
                best = row
                (ROOT / "input_dihedral14_best.json").write_text(
                    json.dumps(candidate, indent=2) + "\n")
    result = {"trials": len(rows),
              "minimum_cnot": min(row["cnot_count"] for row in rows),
              "best": best, "rows": rows}
    (ROOT / "input_dihedral14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "minimum_cnot", "best")}, indent=2))


if __name__ == "__main__":
    main()
