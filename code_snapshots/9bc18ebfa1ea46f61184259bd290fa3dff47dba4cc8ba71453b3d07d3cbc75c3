"""Bounded two-output-CNOT gauge scan of the exact14 circuit.

Two appended CNOTs permute eigenrows and preserve diagonalization. Recompile
all nontrivial ordered pairs with one fixed Qiskit level-3 setting.
"""

import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
EDGES = [(a, b) for a in range(4) for b in range(4) if a != b]


def main():
    original = to_qiskit(json.loads(SOURCE.read_text()))
    rows = []
    best = None
    for edge1 in EDGES:
        for edge2 in EDGES:
            if edge1 == edge2:
                continue
            circuit = original.copy()
            for a, b in (edge1, edge2):
                circuit.cx(3 - a, 3 - b)
            compiled = transpile(circuit, basis_gates=["u", "cx"],
                                 optimization_level=3, seed_transpiler=42)
            candidate = from_qiskit(compiled)
            check = evaluate(candidate)
            row = {"output_cnots": [edge1, edge2],
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid_diagonalizer": check["valid_diagonalizer"]}
            rows.append(row)
            if check["valid_diagonalizer"] and (
                    best is None or row["cnot_count"] < best[0]["cnot_count"]):
                best = row, candidate
    if best:
        (ROOT / "search13_two_output_cx_qiskit_best.json").write_text(
            json.dumps(best[1], indent=2) + "\n")
    print(json.dumps({"source": SOURCE.name, "seed": 42, "rows": rows,
                      "best": best[0] if best else None}, indent=2))


if __name__ == "__main__":
    main()
