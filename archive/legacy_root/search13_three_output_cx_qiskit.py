"""All noncancelling three-CNOT output eigenrow gauges of exact14.

Any output CNOT sequence permutes eigenrows and preserves V4
diagonalization. This is a fixed-seed compiler screen, not an optimum proof.
"""

import json
from collections import Counter
from itertools import product
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
EDGES = [(a, b) for a in range(4) for b in range(4) if a != b]


def main():
    original = to_qiskit(json.loads(SOURCE.read_text()))
    distribution = Counter()
    valid_count = 0
    best = None
    checked = 0
    for sequence in product(EDGES, repeat=3):
        if sequence[0] == sequence[1] or sequence[1] == sequence[2]:
            continue
        circuit = original.copy()
        for a, b in sequence:
            circuit.cx(3 - a, 3 - b)
        compiled = transpile(circuit, basis_gates=["u", "cx"],
                             optimization_level=3, seed_transpiler=42)
        candidate = from_qiskit(compiled)
        check = evaluate(candidate)
        checked += 1
        count = check["cnot_count"]
        distribution[count] += 1
        valid_count += int(check["valid_diagonalizer"])
        if check["valid_diagonalizer"] and (
                best is None or count < best[0]["cnot_count"]):
            best = ({"output_cnots": sequence, "cnot_count": count,
                     "off_diagonal_error": check["off_diagonal_error"]},
                    candidate)
    if best:
        (ROOT / "search13_three_output_cx_qiskit_best.json").write_text(
            json.dumps(best[1], indent=2) + "\n")
    result = {"source": SOURCE.name, "seed": 42, "checked": checked,
              "valid_count": valid_count,
              "cnot_distribution": dict(sorted(distribution.items())),
              "best": best[0] if best else None}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
