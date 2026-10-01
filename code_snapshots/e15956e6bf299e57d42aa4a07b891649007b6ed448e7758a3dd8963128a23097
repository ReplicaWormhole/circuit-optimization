"""Bounded Qiskit recompilation of one-output-CNOT gauges of exact14.

An output CNOT only permutes computational-basis eigenrows, so every input
to this scan remains an exact diagonalizer. Compiler results are numerical
evidence and require an exact follow-up if a shorter circuit appears.
"""

import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
SEED = 42


def main():
    source = json.loads(SOURCE.read_text())
    original = to_qiskit(source)
    rows = []
    best = None
    for control in range(4):
        for target in range(4):
            if control == target:
                continue
            trial = original.copy()
            trial.cx(3 - control, 3 - target)
            compiled = transpile(trial, basis_gates=["u", "cx"],
                                 optimization_level=3,
                                 seed_transpiler=SEED)
            candidate = from_qiskit(compiled)
            check = evaluate(candidate)
            row = {
                "output_cnot": [control, target],
                "compiled_cnot_count": check["cnot_count"],
                "off_diagonal_error": check["off_diagonal_error"],
                "valid_diagonalizer": check["valid_diagonalizer"],
            }
            rows.append(row)
            if check["valid_diagonalizer"] and (
                    best is None or row["compiled_cnot_count"]
                    < best[0]["compiled_cnot_count"]):
                best = row, candidate
    if best:
        (ROOT / "search13_output_cx_qiskit_best.json").write_text(
            json.dumps(best[1], indent=2) + "\n")
    print(json.dumps({"source": SOURCE.name, "seed": SEED,
                      "rows": rows, "best": best[0] if best else None}, indent=2))


if __name__ == "__main__":
    main()
