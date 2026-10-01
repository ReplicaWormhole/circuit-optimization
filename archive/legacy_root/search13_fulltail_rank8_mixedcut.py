"""Mixed input/output tensor-cut screen for ordered seven-CNOT tails.

A circuit tensor network has bond dimension two on every wire segment.
For each boundary bipartition, exhaustive assignment of its seven CNOT
vertices gives a valid mincut rank cap. Compare with the exact-gauge
tail's numerical realignment ranks; exactify a violation separately.
"""

from itertools import product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from search13_fulltail_rank8_fit import SCHEDULES
from search13_fulltail_invariant_gauge_rank import setup

ROOT = Path(__file__).resolve().parent
WITNESS = ROOT / "search13_fulltail_invariant_exact_witness_result.json"
OUT = ROOT / "search13_fulltail_rank8_mixedcut_result.json"


def exact_gauge_numeric():
    data = json.loads(WITNESS.read_text())
    zeta = sp.exp(2 * sp.pi * sp.I / 96)
    g = np.zeros((16, 16), dtype=complex)
    for row, columns in data["gauge_sparse_rows"].items():
        for col, expression in columns.items():
            g[int(row), int(col)] = complex(sp.N(sp.sympify(
                expression, locals={"zeta": zeta}), 17))
    return g


def wire_graph(schedule):
    # Boundary vertices 0..3 are outputs and 4..7 are inputs; CNOT
    # vertices 8..14. The tensor axes use that same boundary ordering.
    edges = []
    for wire in range(4):
        incidents = [8 + t for t, pair in enumerate(schedule) if wire in pair]
        path = [4 + wire] + incidents + [wire]
        edges.extend(zip(path, path[1:]))
    return edges


def mincut(edges, boundary_mask):
    best = 99
    for gate_mask in range(1 << 7):
        sides = [(boundary_mask >> j) & 1 for j in range(8)]
        sides.extend((gate_mask >> j) & 1 for j in range(7))
        cut = sum(sides[a] != sides[b] for a, b in edges)
        best = min(best, cut)
    return best


def realignment_rank(matrix, boundary_mask):
    left = tuple(j for j in range(8) if boundary_mask & (1 << j))
    right = tuple(j for j in range(8) if j not in left)
    a = matrix.reshape((2,) * 8).transpose(left + right)
    s = np.linalg.svd(a.reshape(1 << len(left), 1 << len(right)),
                      compute_uv=False)
    return int(np.count_nonzero(s > 1e-9)), float(s[-1])


def main():
    tail, _, _, _ = setup()
    target = exact_gauge_numeric() @ tail
    ranks = {mask: realignment_rank(target, mask)
             for mask in range(1, 256, 2) if mask != 255}
    rows = []
    for case, schedule in enumerate(SCHEDULES):
        edges = wire_graph(schedule)
        violations = []
        for mask, (rank, least) in ranks.items():
            cut = mincut(edges, mask)
            if rank > 2**cut:
                violations.append({"mask": mask,
                                   "left_boundary_axes": [j for j in range(8)
                                                          if mask & (1 << j)],
                                   "numerical_rank": rank,
                                   "mincut_edges": cut, "rank_cap": 2**cut,
                                   "smallest_singular": least})
        rows.append({"case": case, "schedule": schedule,
                     "tested_boundary_cuts": len(ranks),
                     "violation_count": len(violations),
                     "violations": violations})
        print(json.dumps({"case": case, "violations": len(violations),
                          "first": violations[:1]}), flush=True)
    result = {"exact_gauge_source": WITNESS.name,
              "boundary_axis_order": "output0..3,input0..3",
              "rank_tolerance": 1e-9,
              "rows": rows,
              "scope": "numerical target ranks with exact graph mincuts; exactify any violation before claiming obstruction"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
