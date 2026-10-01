"""Rank-screen three-cycle eigenrow permutations of exact15's final suffix.

A row permutation within one V4 eigenvalue sector commutes with the output
diagonal. Test whether it changes the fixed ten-CNOT-prefix final suffix so
that a four-CNOT replacement is not excluded by operator-Schmidt cut ranks.
This gives necessary conditions only; it is not a synthesis or lower bound
for arbitrary diagonalizers.
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np

from exact_check import exact_eigenvalue_labels
from topology14_alt_suffix_rank import CUTS, EDGES
from topology14_boundary_schmidt import schmidt_singular_values
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"
RESULT = ROOT / "output_threecycle_suffix14_result.json"


def ranks_of(matrix):
    return {"".join(map(str, cut)): int(np.count_nonzero(
        schmidt_singular_values(matrix, cut) > 1e-9)) for cut in CUTS}


def passing_graphs(ranks):
    required = {cut: math.ceil(math.log2(rank))
                for cut, rank in ranks.items()}
    good = []
    for graph in itertools.combinations_with_replacement(EDGES, 4):
        if all(sum((a in CUTS[index]) != (b in CUTS[index])
                   for a, b in graph) >= required["".join(map(str, CUTS[index]))]
               for index in range(len(CUTS))):
            good.append([[a, b] for a, b in graph])
    return good


def main():
    source = json.loads(SOURCE.read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    gates = source["gates"]
    cx_positions = [i for i, gate in enumerate(gates) if gate["gate"] == "cx"]
    assert len(cx_positions) == 15
    split = cx_positions[-5]
    suffix = circuit_matrix(gates[split:])
    rows = []
    for label in sorted(set(labels)):
        sector = [i for i, value in enumerate(labels) if value == label]
        for triple in itertools.combinations(sector, 3):
            a, b, c = triple
            for direction in (0, 1):
                perm = list(range(16))
                perm[a], perm[b], perm[c] = ((b, c, a) if direction == 0
                                               else (c, a, b))
                changed = np.eye(16, dtype=complex)[perm] @ suffix
                ranks = ranks_of(changed)
                graphs = passing_graphs(ranks)
                rows.append({"label": label, "triple": triple,
                             "direction": direction, "ranks": ranks,
                             "compatible_four_cx_graphs": graphs})
    result = {"trials": len(rows), "split_gate": split,
              "baseline_ranks": ranks_of(suffix),
              "compatible_trials": sum(bool(r["compatible_four_cx_graphs"])
                                       for r in rows),
              "rows": rows,
              "rank_tolerance": 1e-9,
              "scope": "single three-cycle row permutation in one eigenvalue sector; fixed ten-CNOT prefix"}
    RESULT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"trials": result["trials"],
                      "baseline_ranks": result["baseline_ranks"],
                      "compatible_trials": result["compatible_trials"]},
                     indent=2))


if __name__ == "__main__":
    main()
