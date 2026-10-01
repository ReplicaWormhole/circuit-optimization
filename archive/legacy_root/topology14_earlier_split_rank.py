"""Find earliest suffix cut where one-CNOT reduction passes rank-cut filters.

For each of four 15-CNOT catalog circuits, isolate the final k CNOTs for
k=5..9. The proposed replacement uses k-1 CNOTs. Enumerate all undirected
multigraphs with k-1 edges and test seven operator-Schmidt cut constraints.
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


def valid_graphs(bounds, edge_count):
    result = []
    for graph in itertools.combinations_with_replacement(EDGES, edge_count):
        if all(sum((a in cut) != (b in cut) for a, b in graph) >= bound
               for cut, bound in bounds):
            result.append(graph)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for index in range(4):
        path = ROOT / f"fivemask15_topology_{index}.json"
        gates = json.loads(path.read_text())["gates"]
        cx_positions = [i for i, gate in enumerate(gates) if gate["gate"] == "cx"]
        for count in range(5, 10):
            start = cx_positions[-count]
            suffix = gates[start:]
            u = circuit_matrix(suffix)
            cuts = []
            for cut in CUTS:
                singular = schmidt_singular_values(u, cut)
                rank = int(np.count_nonzero(singular > 1e-9))
                cuts.append({"cut": "".join(map(str, cut)), "rank": rank,
                             "bound": math.ceil(math.log2(rank)),
                             "smallest_nonzero_singular": float(singular[rank - 1])})
            graphs = valid_graphs([(cut, row["bound"]) for cut, row in zip(CUTS, cuts)],
                                  count - 1)
            rows.append({"candidate": path.name, "old_suffix_cnot": count,
                         "new_suffix_cnot": count - 1, "start_gate": start,
                         "original_suffix_cnot_edges": [[g["control"], g["target"]]
                                                       for g in suffix if g["gate"] == "cx"],
                         "cuts": cuts, "passing_graph_count": len(graphs),
                         "example_graphs": [[list(edge) for edge in graph]
                                            for graph in graphs[:10]]})
    result = {"results": rows, "necessary_only": True, "rank_tolerance": 1e-9}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"results": [{"candidate": row["candidate"],
                                    "old_suffix_cnot": row["old_suffix_cnot"],
                                    "passing_graph_count": row["passing_graph_count"],
                                    "ranks": {cut["cut"]: cut["rank"]
                                              for cut in row["cuts"]}}
                                   for row in rows]}, indent=2))


if __name__ == "__main__":
    main()
