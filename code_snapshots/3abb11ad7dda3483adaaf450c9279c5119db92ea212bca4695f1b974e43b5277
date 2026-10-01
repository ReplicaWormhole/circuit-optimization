"""Screen products of disjoint nonlocal output Givens gauges.

Each factor mixes two computational output rows carrying the same V4
eigenvalue, so their disjoint product commutes with the output eigenvalue
matrix.  We test a finite exact-angle grid for a four-CNOT replacement of
the final five-CNOT suffix after the exact15 ten-CNOT prefix.
"""

import argparse
import collections
import itertools
import json
import math
from pathlib import Path

from exact_check import exact_eigenvalue_labels
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block
from delete14_nonlocal_gauge_screen import (CUTS, feasible_graphs,
                                           schmidt_rank, two_level, unitary)


ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    split = 13 + len(exact_block()) + len(reverse_block())
    suffix = unitary(source["gates"][split:])
    pairs = [(a, b) for a, b in itertools.combinations(range(16), 2)
             if labels[a] == labels[b] and (a ^ b).bit_count() >= 2]
    records = []
    compatible = []
    signatures = collections.Counter()
    phases = [k * math.pi / 2 for k in range(4)]
    for first, second in itertools.combinations(pairs, 2):
        if set(first) & set(second):
            continue
        for first_phase_index, first_phase in enumerate(phases):
            first_gauge = two_level(*first, math.pi / 2, first_phase)
            for second_phase_index, second_phase in enumerate(phases):
                gauge = (two_level(*second, math.pi / 2, second_phase)
                         @ first_gauge)
                ranks = [schmidt_rank(gauge @ suffix, cut)[0] for cut in CUTS]
                bounds = [int(math.ceil(math.log2(rank))) for rank in ranks]
                graphs = feasible_graphs(bounds)
                row = {"first_pair": list(first), "second_pair": list(second),
                       "first_phase_over_pi": first_phase_index / 2,
                       "second_phase_over_pi": second_phase_index / 2,
                       "angles_over_pi": [0.5, 0.5],
                       "ranks": ranks, "cut_cnot_bounds": bounds,
                       "feasible_graph_count": len(graphs)}
                records.append(row)
                signatures[str(ranks)] += 1
                if graphs:
                    compatible.append({**row, "feasible_graphs": graphs})
    result = {"source": "topology16_15_exact_candidate.json",
              "split_gate": split, "prefix_cnot": 10,
              "pairs_available": len(pairs),
              "disjoint_pair_choices": len(records) // 16,
              "trials": len(records), "compatible_count": len(compatible),
              "rank_signature_histogram": dict(signatures),
              "compatible": compatible, "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("records", "compatible")}, indent=2))
    print(json.dumps({"first_compatible": compatible[:5]}, indent=2))


if __name__ == "__main__":
    main()
