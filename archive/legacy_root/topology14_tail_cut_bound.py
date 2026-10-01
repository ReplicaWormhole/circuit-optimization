"""Rank-cut constraints for a four-CNOT realization of the fixed final suffix.

For each four-qubit bipartition, operator-Schmidt rank r requires at least
ceil(log2 r) crossing CNOTs. Enumerate all four-edge undirected multigraphs
and retain those satisfying all seven cut constraints. This is a necessary
topology filter, not a synthesis test or an absolute lower bound.
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
EDGES = tuple(itertools.combinations(range(4), 2))
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    suffix = baseline["gates"][38:]
    u = circuit_matrix(suffix)
    ranks = {}
    for cut in CUTS:
        singular = schmidt_singular_values(u, cut)
        rank = int(np.count_nonzero(singular > 1e-9))
        ranks["".join(map(str, cut))] = {"rank": rank,
                                        "min_crossing_cnot": math.ceil(math.log2(rank)),
                                        "smallest_nonzero_singular": float(singular[rank - 1])}
    passing = []
    for graph in itertools.combinations_with_replacement(EDGES, 4):
        if all(sum((edge[0] in cut) != (edge[1] in cut) for edge in graph)
               >= ranks["".join(map(str, cut))]["min_crossing_cnot"]
               for cut in CUTS):
            passing.append(graph)
    result = {"original_directed_cnot_edges": [[g["control"], g["target"]]
                                               for g in suffix if g["gate"] == "cx"],
              "cut_ranks": ranks, "four_edge_graphs_total": 126,
              "passing_four_edge_graphs": [[list(edge) for edge in graph]
                                            for graph in passing],
              "passing_count": len(passing)}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "passing_four_edge_graphs"}, indent=2))


if __name__ == "__main__":
    main()
