"""Rank-screen simultaneous three-cycle output gauges in two eigensectors.

Run 139 rejected every single three-cycle row permutation for a fixed
ten-CNOT prefix. This changes multiple eigensectors coherently by composing
two disjoint three-cycles with distinct eigenvalue labels. It retains the
same fixed prefix and only tests necessary operator-Schmidt cut bounds.
"""

import itertools
import json
from pathlib import Path

import numpy as np

from exact_check import exact_eigenvalue_labels
from output_threecycle_suffix14 import ranks_of, passing_graphs
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"
RESULT = ROOT / "output_two_threecycles_suffix14_result.json"


def cycles(labels):
    groups = {}
    for label in sorted(set(labels)):
        sector = [i for i, value in enumerate(labels) if value == label]
        groups[label] = []
        for triple in itertools.combinations(sector, 3):
            a, b, c = triple
            for image in ((b, c, a), (c, a, b)):
                groups[label].append((triple, image))
    return groups


def main():
    source = json.loads(SOURCE.read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    gates = source["gates"]
    cx_positions = [i for i, gate in enumerate(gates) if gate["gate"] == "cx"]
    assert len(cx_positions) == 15
    suffix = circuit_matrix(gates[cx_positions[-5]:])
    groups = cycles(labels)
    rows = []
    for left, right in itertools.combinations(groups, 2):
        for (triple1, image1), (triple2, image2) in itertools.product(
                groups[left], groups[right]):
            perm = list(range(16))
            for triple, image in ((triple1, image1), (triple2, image2)):
                for source_row, target_row in zip(triple, image):
                    perm[source_row] = target_row
            changed = np.eye(16, dtype=complex)[perm] @ suffix
            ranks = ranks_of(changed)
            graphs = passing_graphs(ranks)
            rows.append({"labels": [left, right],
                         "triples": [triple1, triple2],
                         "images": [image1, image2], "ranks": ranks,
                         "compatible_four_cx_graphs": graphs})
    result = {"trials": len(rows),
              "compatible_trials": sum(bool(row["compatible_four_cx_graphs"])
                                       for row in rows),
              "best_crosscut_ranks": {
                  cut: min(row["ranks"][cut] for row in rows)
                  for cut in ("01", "02", "03")},
              "rows": rows, "rank_tolerance": 1e-9,
              "scope": "two three-cycles in distinct eigenvalue sectors; fixed ten-CNOT prefix"}
    RESULT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "trials", "compatible_trials", "best_crosscut_ranks")}, indent=2))


if __name__ == "__main__":
    main()
