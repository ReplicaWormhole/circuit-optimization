"""Scan all output wire permutations and bit flips of the 15-CNOT circuit.

These 384 monomial output maps preserve diagonalization while reordering the
measurement outcomes. Qiskit resynthesizes the complete circuit. This does
not cover arbitrary affine maps or rotations within degenerate eigenspaces.
"""

import itertools
import json
from pathlib import Path

from qiskit import transpile

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
    for permutation in itertools.permutations(range(4)):
        swaps = swaps_for(permutation)
        for flips in range(16):
            circuit = to_qiskit(base)
            for a, b in swaps:
                circuit.swap(3 - a, 3 - b)
            for wire in range(4):
                if flips & (1 << wire):
                    circuit.x(3 - wire)
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = from_qiskit(optimized)
            check = evaluate(candidate)
            row = {"permutation": permutation, "flips": flips,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid": check["valid_diagonalizer"]}
            rows.append(row)
            if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                    best["cnot_count"], best["off_diagonal_error"]):
                best = row
                (ROOT / "output_wire_permutation14_best.json").write_text(
                    json.dumps(candidate, indent=2) + "\n")
    result = {"trials": len(rows),
              "minimum_cnot": min(row["cnot_count"] for row in rows),
              "best": best, "rows": rows}
    (ROOT / "output_wire_permutation14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "minimum_cnot", "best")}, indent=2))


if __name__ == "__main__":
    main()
