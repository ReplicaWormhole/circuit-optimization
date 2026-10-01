"""Enumerate binary circulant input symmetries and synthesize non-shift cases.

For A=p(V4) over F2, A commutes exactly with the cyclic bit shift V4.  If A
is invertible, precomposing a cycle diagonalizer with the corresponding CNOT
permutation preserves its diagonalization.  BFS finds shortest CNOT words for
the four non-shift units; Qiskit gives bounded synthesis upper bounds.
"""

import argparse
from collections import deque
import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
IDENTITY = (8, 4, 2, 1)
SHIFT = (1, 8, 4, 2)
EDGES = [(c, t) for c in range(4) for t in range(4) if c != t]


def compose(a, b):
    return tuple(xor_rows(b, mask) for mask in a)


def xor_rows(rows, mask):
    value = 0
    for q in range(4):
        if mask & (1 << (3 - q)):
            value ^= rows[q]
    return value


def polynomial_map(poly):
    powers = [IDENTITY]
    for _ in range(3):
        powers.append(compose(SHIFT, powers[-1]))
    return tuple(_polynomial_row(poly, powers, q) for q in range(4))


def _polynomial_row(poly, powers, q):
    value = 0
    for k in range(4):
        if poly & (1 << k):
            value ^= powers[k][q]
    return value


def apply_cx(rows, edge):
    c, t = edge
    out = list(rows)
    out[t] ^= out[c]
    return tuple(out)


def shortest_words(targets):
    todo = deque([IDENTITY])
    parent = {IDENTITY: None}
    wanted = set(targets)
    while todo and wanted:
        current = todo.popleft()
        wanted.discard(current)
        for edge in EDGES:
            nxt = apply_cx(current, edge)
            if nxt not in parent:
                parent[nxt] = (current, edge)
                todo.append(nxt)
    if wanted:
        raise RuntimeError(f"unreachable invertible maps: {wanted}")
    words = {}
    for target in targets:
        word = []
        current = target
        while current != IDENTITY:
            previous, edge = parent[current]
            word.append(edge)
            current = previous
        words[target] = list(reversed(word))
    return words


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    args = p.parse_args()
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    unit_polys = [poly for poly in range(16) if poly.bit_count() % 2 == 1]
    maps = {poly: polynomial_map(poly) for poly in unit_polys}
    assert len(set(maps.values())) == 8
    assert all(compose(matrix, SHIFT) == compose(SHIFT, matrix)
               for matrix in maps.values())
    words = shortest_words(maps.values())
    shifts = {maps[1 << k] for k in range(4)}
    results = []
    best = None
    best_candidate = None
    for poly, matrix in maps.items():
        word = words[matrix]
        row = {"polynomial_bits": poly, "matrix_rows": matrix,
               "minimum_input_cnot": len(word), "input_cx": word,
               "is_shift_power": matrix in shifts}
        if matrix not in shifts:
            gates = [{"gate": "cx", "control": c, "target": t}
                     for c, t in word]
            circuit = to_qiskit({"n": 4, "gates": gates + base["gates"]})
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3,
                                  seed_transpiler=args.seed)
            candidate = from_qiskit(optimized)
            check = evaluate(candidate)
            row.update(total_cnot_count=check["cnot_count"],
                       off_diagonal_error=check["off_diagonal_error"],
                       valid=check["valid_diagonalizer"])
            if best is None or (row["total_cnot_count"],
                                row["off_diagonal_error"]) < (
                                    best["total_cnot_count"],
                                    best["off_diagonal_error"]):
                best, best_candidate = row, candidate
        results.append(row)
    args.output.write_text(json.dumps(best_candidate, indent=2) + "\n")
    summary = {"exact_commuting_maps": len(maps),
               "non_shift_maps_synthesized": len(maps) - len(shifts),
               "seed": args.seed, "best": best, "results": results}
    args.result.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
