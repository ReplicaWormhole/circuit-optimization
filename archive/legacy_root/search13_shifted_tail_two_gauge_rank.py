"""Exact modular cut-rank screen for gauged full eight-CNOT matchgate tail.

Shift the synthesis boundary immediately after the exact six-CNOT prefix.
The fixed tail now includes both first-layer pair blocks, the center local
layer, and both terminal pair blocks. For pairs of cross-terminal masks and
angles k*pi/16 with k in {-4,-2,2,4}, compute exact lower cut ranks in
F193 and F769 and enumerate minimal cut-compatible CNOT multigraphs.
Passing a seven-CNOT graph screen is necessary, never sufficient.
"""

import itertools
import json
from pathlib import Path

import numpy as np

from algebraic13_crosspair_parity_gauge_rank import (
    CUTS, EDGES, MASKS, minimum_graph_edges, rank_mod, realign,
)
from delete14_integer_audit import matrix


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
PRIMES = (193, 769)
ANGLES = (-4, -2, 2, 4)


def primitive_96_root(p):
    for z in range(2, p):
        if pow(z, 96, p) == 1 and pow(z, 48, p) != 1 and pow(z, 32, p) != 1:
            return z
    raise AssertionError(f"no primitive 96th root in F_{p}")


def gate_mod(gate, p, powers):
    if gate["gate"] == "cx":
        result = np.zeros((16, 16), dtype=np.int64)
        c, t = gate["control"], gate["target"]
        for col in range(16):
            row = col ^ ((8 >> t) if col & (8 >> c) else 0)
            result[row, col] = 1
        return result
    q = gate["qubit"]
    local, _ = matrix(gate)
    result = np.zeros((16, 16), dtype=np.int64)
    for col in range(16):
        b = int(bool(col & (8 >> q)))
        for a in (0, 1):
            row = (col & ~(8 >> q)) | (a * (8 >> q))
            result[row, col] = int(np.dot(powers, np.asarray(local[a][b],
                                                             dtype=np.int64)) % p)
    return result


def tail_mod(gates, p, z):
    powers = np.array([pow(z, k, p) for k in range(32)], dtype=np.int64)
    result = np.eye(16, dtype=np.int64)
    for gate in gates:
        result = gate_mod(gate, p, powers) @ result % p
    return result


def parity_sign(basis, mask):
    return 1 if (basis & mask).bit_count() % 2 == 0 else -1


def main():
    gates = json.loads(BASE.read_text())["gates"]
    tail = gates[13:]
    assert sum(g["gate"] == "cx" for g in gates[:13]) == 6
    assert sum(g["gate"] == "cx" for g in tail) == 8
    fields = [(p, primitive_96_root(p)) for p in PRIMES]
    tails = {p: tail_mod(tail, p, z) for p, z in fields}
    rows = []
    for first, second in itertools.combinations(MASKS, 2):
        for k1, k2 in itertools.product(ANGLES, repeat=2):
            rank_arrays = []
            for p, z in fields:
                z32 = pow(z, 3, p)
                phases = np.array([pow(z32,
                                       k1 * parity_sign(basis, first)
                                       + k2 * parity_sign(basis, second), p)
                                   for basis in range(16)], dtype=np.int64)
                u = phases[:, None] * tails[p] % p
                rank_arrays.append([rank_mod(realign(u, cut), p)
                                    for cut in CUTS])
            certified = [max(ranks) for ranks in zip(*rank_arrays)]
            bits = [(rank - 1).bit_length() for rank in certified]
            graph_min, witness = minimum_graph_edges(np.asarray(bits))
            rows.append({"masks": [first, second], "angle_numerators": [k1, k2],
                         "modular_ranks_by_field": rank_arrays,
                         "certified_rank_lower_bounds": certified,
                         "cut_cnot_lower_bounds": bits,
                         "minimum_graph_edges": graph_min,
                         "feasible_edge_multiplicities": witness})
    result = {"boundary_gate_start": 13, "base": BASE.name,
              "tail_cnot_count": 8, "fields": [{"prime": p,
              "primitive_96_root": z} for p, z in fields],
              "cuts": [list(c) for c in CUTS],
              "edge_order": [list(e) for e in EDGES],
              "masks": list(MASKS), "angle_numerators": list(ANGLES),
              "rows": rows}
    path = ROOT / "search13_shifted_tail_two_gauge_rank_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "cases": len(rows),
                      "graph_min_histogram": {str(n): sum(
                          r["minimum_graph_edges"] == n for r in rows)
                          for n in range(9)},
                      "seven_cnot_compatible": sum(
                          r["minimum_graph_edges"] <= 7 for r in rows),
                      "best_graph_lower_bound": min(
                          r["minimum_graph_edges"] for r in rows)}, indent=2))


if __name__ == "__main__":
    main()
