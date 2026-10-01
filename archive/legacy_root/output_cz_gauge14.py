"""Search output-diagonal CZ gauges of the exact 15-CNOT V4 circuit.

Every product of output CZ gates commutes with the measured diagonal D, so
it preserves diagonalization while changing the synthesized eigenbasis.
This is a bounded 64-point Clifford gauge scan; it does not search arbitrary
output-local or degenerate-eigenspace rotations.
"""

import itertools
import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
PAIRS = tuple(itertools.combinations(range(4), 2))


def cz_gates(pair):
    control, target = pair
    return [{"gate": "h", "qubit": target},
            {"gate": "cx", "control": control, "target": target},
            {"gate": "h", "qubit": target}]


def main():
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    results = []
    best = None
    for mask in range(1 << len(PAIRS)):
        extra = [gate for index, pair in enumerate(PAIRS) if mask & (1 << index)
                 for gate in cz_gates(pair)]
        input_candidate = {"n": 4, "gates": base["gates"] + extra}
        optimized = transpile(to_qiskit(input_candidate),
                              basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        candidate = from_qiskit(optimized)
        check = evaluate(candidate)
        row = {"mask": mask, "pairs": [pair for index, pair in enumerate(PAIRS)
                                         if mask & (1 << index)],
               "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid": check["valid_diagonalizer"]}
        results.append(row)
        if best is None or (row["cnot_count"], row["off_diagonal_error"]) < (
                best["cnot_count"], best["off_diagonal_error"]):
            best = row
            (ROOT / "output_cz_gauge14_best.json").write_text(
                json.dumps(candidate, indent=2) + "\n")
    summary = {"pairs": PAIRS, "gauges_tested": len(results),
               "minimum_cnot": min(row["cnot_count"] for row in results),
               "valid_count": sum(row["valid"] for row in results),
               "best": best, "results": results}
    (ROOT / "output_cz_gauge14_result.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in (
        "gauges_tested", "minimum_cnot", "valid_count", "best")}, indent=2))


if __name__ == "__main__":
    main()
