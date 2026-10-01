"""Recompile exact14 after every binary-linear input symmetry of V4.

The eight invertible 4x4 circulant binary maps commute with the right shift.
Shortest CNOT implementations are found by BFS over GL(4,2). Each map is
prepended to the exact14 diagonalizer and recompiled by fixed-seed Qiskit.
"""

import json
from collections import deque
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
IDENTITY = (8, 4, 2, 1)
EDGES = [(a, b) for a in range(4) for b in range(4) if a != b]


def circulant(coefficients):
    return tuple(sum(1 << (3 - j) for j in range(4)
                     if (coefficients >> ((i - j) % 4)) & 1)
                 for i in range(4))


def shortest_paths(targets):
    queue = deque([IDENTITY])
    paths = {IDENTITY: ()}
    while queue and not targets.issubset(paths):
        masks = queue.popleft()
        for control, target in EDGES:
            nxt = list(masks)
            nxt[target] ^= nxt[control]
            nxt = tuple(nxt)
            if nxt not in paths:
                paths[nxt] = paths[masks] + ((control, target),)
                queue.append(nxt)
    assert targets.issubset(paths)
    return paths


def main():
    source = json.loads(SOURCE.read_text())
    targets = {circulant(a) for a in range(16) if a.bit_count() % 2 == 1}
    assert len(targets) == 8 and IDENTITY in targets
    paths = shortest_paths(targets)
    rows = []
    best = None
    for coefficients in range(16):
        if coefficients.bit_count() % 2 != 1:
            continue
        linear = circulant(coefficients)
        prefix = [{"gate": "cx", "control": c, "target": t}
                  for c, t in paths[linear]]
        candidate = {"n": 4, "gates": prefix + source["gates"]}
        circuit = to_qiskit(candidate)
        compiled = transpile(circuit, basis_gates=["u", "cx"],
                             optimization_level=3, seed_transpiler=42)
        result = from_qiskit(compiled)
        check = evaluate(result)
        row = {"coefficients": coefficients, "matrix_rows": linear,
               "input_cnot_count": len(prefix),
               "compiled_cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        if check["valid_diagonalizer"] and (
                best is None or row["compiled_cnot_count"]
                < best[0]["compiled_cnot_count"]):
            best = row, result
    if best:
        (ROOT / "search13_circulant_input_gauge_best.json").write_text(
            json.dumps(best[1], indent=2) + "\n")
    print(json.dumps({"source": SOURCE.name, "seed": 42,
                      "rows": rows, "best": best[0] if best else None}, indent=2))


if __name__ == "__main__":
    main()
