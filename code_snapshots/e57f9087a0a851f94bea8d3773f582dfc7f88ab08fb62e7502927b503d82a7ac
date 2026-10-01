"""Bounded Qiskit pass-order/basis scan of the exact 15-CNOT circuit.

This tests whether transpiler seed and one-qubit basis choices expose a
two-qubit cancellation missed by the seed-42 U/CX compilation. It is a
heuristic compiler search, not an exact synthesis theorem.
"""

import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
BASES = (("u", "cx"), ("rz", "sx", "x", "cx"))
SEEDS = range(64)


def main():
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    circuit = to_qiskit(source)
    rows = []
    best = None
    for basis in BASES:
        for seed in SEEDS:
            optimized = transpile(circuit, basis_gates=list(basis),
                                  optimization_level=3, seed_transpiler=seed)
            # Normalize to the adapter's supported U/CX gate set.
            normalized = transpile(optimized, basis_gates=["u", "cx"],
                                   optimization_level=0)
            candidate = from_qiskit(normalized)
            check = evaluate(candidate)
            row = {"basis": basis, "seed": seed,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid": check["valid_diagonalizer"]}
            rows.append(row)
            if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                    best["cnot_count"], best["off_diagonal_error"]):
                best = row
                (ROOT / "qiskit_seed_basis14_best.json").write_text(
                    json.dumps(candidate, indent=2) + "\n")
    result = {"trials": len(rows),
              "minimum_cnot": min(row["cnot_count"] for row in rows),
              "best": best, "rows": rows}
    (ROOT / "qiskit_seed_basis14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "minimum_cnot", "best")}, indent=2))


if __name__ == "__main__":
    main()
