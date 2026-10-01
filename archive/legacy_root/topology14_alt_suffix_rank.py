"""Cut-rank topology filter for four alternative 15-CNOT diagonalizers.

At each candidate, isolate the final five-CNOT suffix. Enumerate four-edge
undirected CNOT graphs satisfying its seven operator-Schmidt cut ranks.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from topology14_boundary_schmidt import schmidt_singular_values
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
EDGES = tuple(itertools.combinations(range(4), 2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for index in range(4):
        path = ROOT / f"fivemask15_topology_{index}.json"
        gates = json.loads(path.read_text())["gates"]
        cx_positions = [i for i, gate in enumerate(gates) if gate["gate"] == "cx"]
        assert len(cx_positions) == 15
        suffix = gates[cx_positions[-5]:]
        u = circuit_matrix(suffix)
        ranks = {}
        for cut in CUTS:
            singular = schmidt_singular_values(u, cut)
            rank = int(np.count_nonzero(singular > 1e-9))
            ranks["".join(map(str, cut))] = {
                "rank": rank, "min_crossing_cnot": math.ceil(math.log2(rank)),
                "smallest_nonzero_singular": float(singular[rank - 1])}
        passing = []
        for graph in itertools.combinations_with_replacement(EDGES, 4):
            if all(sum((a in cut) != (b in cut) for a, b in graph)
                   >= ranks["".join(map(str, cut))]["min_crossing_cnot"]
                   for cut in CUTS):
                passing.append(graph)
        rows.append({"candidate": path.name,
                     "suffix_start_gate": cx_positions[-5],
                     "suffix_cnot_edges": [[gate["control"], gate["target"]]
                                           for gate in suffix if gate["gate"] == "cx"],
                     "cut_ranks": ranks,
                     "passing_four_edge_graphs": [[list(edge) for edge in graph]
                                                  for graph in passing],
                     "passing_count": len(passing)})
    result = {"results": rows, "rank_tolerance": 1e-9,
              "scope": "fixed suffix targets; local gates may vary in a synthesis"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"results": [{"candidate": r["candidate"],
                                    "suffix_cnot_edges": r["suffix_cnot_edges"],
                                    "ranks": {k: v["rank"] for k, v in r["cut_ranks"].items()},
                                    "passing_count": r["passing_count"]}
                                   for r in rows]}, indent=2))


if __name__ == "__main__":
    main()
