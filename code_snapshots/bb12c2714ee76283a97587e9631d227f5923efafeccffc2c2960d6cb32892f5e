"""Finite exact-modular screen of two cross-pair output parity phases.

For each unordered pair among the nine cross-terminal-pair parity masks,
combine two phases exp(i*k*pi/16*Z_mask) with k in {-4,-2,-1,1,2,4}.
Compute operator-Schmidt ranks in F193 and F257 for all seven wire cuts,
then the smallest cut-compatible undirected CNOT multigraph. This screen
is restricted to its finite angle grid and the fixed terminal suffix S.
"""

import itertools
import json
from pathlib import Path

import numpy as np

from algebraic13_crosspair_parity_gauge_rank import (
    CUTS, EDGES, MASKS, PRIMES, minimum_graph_edges, primitive_32_root,
    rank_mod, realign, suffix_mod,
)


ROOT = Path(__file__).resolve().parent
ANGLES = (-4, -2, -1, 1, 2, 4)


def parity_sign(basis, mask):
    return 1 if (basis & mask).bit_count() % 2 == 0 else -1


def run():
    fields = [(p, primitive_32_root(p)) for p in PRIMES]
    suffixes = {p: suffix_mod(p, root) for p, root in fields}
    rows = []
    for first, second in itertools.combinations(MASKS, 2):
        for k1, k2 in itertools.product(ANGLES, repeat=2):
            rank_arrays = []
            for p, root in fields:
                phases = np.array([pow(root, k1 * parity_sign(basis, first)
                                       + k2 * parity_sign(basis, second), p)
                                   for basis in range(16)], dtype=np.int64)
                u = phases[:, None] * suffixes[p] % p
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
    result = {"suffix": "F01 tensor F23, F=exp(i*pi*(XX+YY)/8)",
              "gauge": "exp(i*k1*pi/16*Z_m1) exp(i*k2*pi/16*Z_m2)",
              "fields": [{"prime": p, "primitive_32_root": root}
                         for p, root in fields],
              "cuts": [list(c) for c in CUTS],
              "edge_order": [list(e) for e in EDGES],
              "masks": list(MASKS), "angle_numerators": list(ANGLES),
              "rows": rows,
              "three_cnot_compatible": [
                  {"masks": r["masks"], "angle_numerators": r["angle_numerators"],
                   "ranks": r["certified_rank_lower_bounds"]}
                  for r in rows if r["minimum_graph_edges"] <= 3]}
    path = ROOT / "algebraic13_crosspair_two_mask_grid_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "cases": len(rows),
                      "three_cnot_compatible": result["three_cnot_compatible"],
                      "graph_min_histogram": {str(n): sum(
                          r["minimum_graph_edges"] == n for r in rows)
                          for n in range(9)},
                      "balanced_cut_rank_minima": [min(
                          r["certified_rank_lower_bounds"][j] for r in rows)
                          for j in (5, 6)]}, indent=2))


if __name__ == "__main__":
    run()
