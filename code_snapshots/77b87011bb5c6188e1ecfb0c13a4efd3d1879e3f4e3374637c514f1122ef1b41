"""Screen one nonlocal output gauge on an alternative nine-CNOT prefix.

The five-mask topology_2 circuit gives a numerical 15-CNOT diagonalizer.
Freeze its first nine CNOTs and ask whether an output SU2 rotation between
equal-eigenvalue rows makes a five-CNOT replacement of its six-CNOT suffix
compatible with all seven operator-Schmidt cut bounds.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import right_shift
from delete14_nonlocal_gauge_screen import (CUTS, PAIRS, schmidt_rank,
                                           two_level, unitary)


ROOT = Path(__file__).resolve().parent
GRAPHS = list(itertools.combinations_with_replacement(PAIRS, 5))


def labels(u):
    diagonal = np.diag(u @ right_shift(4) @ u.conj().T)
    roots = np.array([1, 1j, -1, -1j])
    distance = abs(diagonal[:, None] - roots[None, :])
    if np.max(np.min(distance, axis=1)) > 1e-8:
        raise ValueError("source circuit does not numerically diagonalize V4")
    return np.argmin(distance, axis=1).tolist()


def feasible_graphs(bounds):
    return [graph for graph in GRAPHS
            if all(sum((a in cut) != (b in cut) for a, b in graph) >= bound
                   for cut, bound in zip(CUTS, bounds))]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    source = json.loads(args.source.read_text())
    split = 27
    prefix, suffix = source["gates"][:split], source["gates"][split:]
    if (sum(gate["gate"] == "cx" for gate in prefix) != 9 or
            sum(gate["gate"] == "cx" for gate in suffix) != 6):
        raise ValueError("expected nine-CNOT prefix and six-CNOT suffix")
    output_labels = labels(unitary(source["gates"]))
    suffix_u = unitary(suffix)
    baseline_ranks = [schmidt_rank(suffix_u, cut)[0] for cut in CUTS]
    baseline_bounds = [int(math.ceil(math.log2(rank)))
                       for rank in baseline_ranks]
    baseline_graphs = feasible_graphs(baseline_bounds)
    records = []
    for a, b in itertools.combinations(range(16), 2):
        if output_labels[a] != output_labels[b] or (a ^ b).bit_count() < 2:
            continue
        for angle_num in (1, 2):
            for phase_num in range(4):
                gauge = two_level(a, b, angle_num * math.pi / 4,
                                  phase_num * math.pi / 2)
                ranks = [schmidt_rank(gauge @ suffix_u, cut)[0]
                         for cut in CUTS]
                bounds = [int(math.ceil(math.log2(rank))) for rank in ranks]
                graphs = feasible_graphs(bounds)
                records.append({"row_pair": [a, b],
                                "eigenvalue_label": output_labels[a],
                                "hamming_distance": (a ^ b).bit_count(),
                                "angle_over_pi": angle_num / 4,
                                "phase_over_pi": phase_num / 2,
                                "ranks": ranks, "cut_cnot_bounds": bounds,
                                "feasible_graph_count": len(graphs),
                                "example_graphs": graphs[:16]})
    result = {"source": str(args.source), "split_gate": split,
              "prefix_cnot": 9, "suffix_cnot": 6,
              "cut_order": [list(cut) for cut in CUTS],
              "baseline_ranks": baseline_ranks,
              "baseline_feasible_graph_count": len(baseline_graphs),
              "baseline_graphs": baseline_graphs,
              "gauge_trials": len(records),
              "rank_compatible_trials": sum(row["feasible_graph_count"] > 0
                                            for row in records),
              "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("baseline_graphs", "records")}, indent=2))
    print(json.dumps({"top_graph_counts": sorted(records,
                                                  key=lambda row: -row["feasible_graph_count"])[:5]},
                     indent=2))


if __name__ == "__main__":
    main()
