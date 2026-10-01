"""Scan disjoint two-pair output-eigenspace gauges for boundary rank reduction."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from exact_check import exact_eigenvalue_labels
from topology14_boundary_schmidt import matrix as boundary_matrix
from topology14_boundary_schmidt import schmidt_singular_values
from topology14_twolevel_boundary_rank import circuit_matrix, two_level


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    labels = exact_eigenvalue_labels(baseline)["output_labels"]
    suffix = circuit_matrix(baseline["gates"][38:])
    boundary = boundary_matrix(0.0)
    pairs = [(a, b) for a in range(16) for b in range(a + 1, 16)
             if labels[a] == labels[b]]
    pair_gauges = {}
    for pair in pairs:
        for axis in ("x", "y"):
            for numerator in (1, 2):
                pair_gauges[(pair, axis, numerator)] = two_level(
                    *pair, numerator * math.pi / 4, axis)
    records = []
    for index, first_pair in enumerate(pairs):
        for second_pair in pairs[index + 1:]:
            if set(first_pair) & set(second_pair):
                continue
            for first_axis in ("x", "y"):
                for second_axis in ("x", "y"):
                    for first_num in (1, 2):
                        for second_num in (1, 2):
                            first = pair_gauges[(first_pair, first_axis, first_num)]
                            second = pair_gauges[(second_pair, second_axis, second_num)]
                            changed = suffix.conj().T @ (second @ first) @ suffix @ boundary
                            singular = schmidt_singular_values(changed, (0, 1))
                            records.append({
                                "pairs": [first_pair, second_pair],
                                "axes": [first_axis, second_axis],
                                "angles_over_pi": [f"{first_num}/4", f"{second_num}/4"],
                                "rank_1e-9": int(np.count_nonzero(singular > 1e-9)),
                                "smallest_singular": float(singular[-1]),
                            })
    histogram = {str(rank): sum(row["rank_1e-9"] == rank for row in records)
                 for rank in sorted({row["rank_1e-9"] for row in records})}
    result = {"trial_count": len(records), "rank_histogram": histogram,
              "min_rank": min(row["rank_1e-9"] for row in records),
              "best": sorted(records, key=lambda row: (row["rank_1e-9"],
                                                          row["smallest_singular"]))[:12],
              "scope": "two disjoint equal-label output pairs, X/Y axes, pi/4 or pi/2 angles"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
