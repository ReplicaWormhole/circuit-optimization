"""Scan short output CNOT permutations of the exact 15-CNOT eigenbasis.

Any output computational-basis permutation preserves diagonalization, though
it may permute the reported eigenvalue labels. This bounded scan appends one
or two CNOTs, then asks Qiskit to re-optimize the complete circuit.
"""

import itertools
import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
EDGES = tuple((a, b) for a in range(4) for b in range(4) if a != b)


def main():
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    rows = []
    best = None
    for length in (1, 2):
        for edges in itertools.product(EDGES, repeat=length):
            if length == 2 and edges[0] == edges[1]:
                continue
            extension = [{"gate": "cx", "control": a, "target": b}
                         for a, b in edges]
            circuit = to_qiskit({"n": 4, "gates": base["gates"] + extension})
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = from_qiskit(optimized)
            check = evaluate(candidate)
            row = {"output_edges": edges,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid": check["valid_diagonalizer"]}
            rows.append(row)
            if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                    best["cnot_count"], best["off_diagonal_error"]):
                best = row
                (ROOT / "output_linear_gauge14_best.json").write_text(
                    json.dumps(candidate, indent=2) + "\n")
    result = {"trials": len(rows),
              "minimum_cnot": min(row["cnot_count"] for row in rows),
              "best": best, "rows": rows}
    (ROOT / "output_linear_gauge14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "minimum_cnot", "best")}, indent=2))


if __name__ == "__main__":
    main()
