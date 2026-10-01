"""Screen nonlocal exact output-eigenspace gauges for four-CNOT tails.

An output two-level rotation between rows with the same V4 eigenvalue
commutes with the diagonal eigenvalue matrix.  We multiply it into the final
five-CNOT suffix of the exact15 circuit and test whether any four-edge CNOT
multiset satisfies all seven operator-Schmidt cut lower bounds.  These bounds
are necessary, never sufficient, for a four-CNOT synthesis.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix
from exact_check import exact_eigenvalue_labels
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent
PAIRS = list(itertools.combinations(range(4), 2))
CUTS = [(0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3)]
GRAPHS = list(itertools.combinations_with_replacement(PAIRS, 4))


def unitary(gates):
    result = np.eye(16, dtype=complex)
    for gate in gates:
        matrix = (cnot(gate["control"], gate["target"], 4)
                  if gate["gate"] == "cx" else
                  embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        result = matrix @ result
    return result


def two_level(a, b, angle, phase):
    c, s = math.cos(angle / 2), math.sin(angle / 2)
    result = np.eye(16, dtype=complex)
    result[a, a] = result[b, b] = c
    result[a, b] = -np.exp(-1j * phase) * s
    result[b, a] = np.exp(1j * phase) * s
    return result


def schmidt_rank(matrix, cut):
    other = [q for q in range(4) if q not in cut]
    axes = list(cut) + [q + 4 for q in cut] + other + [q + 4 for q in other]
    shape = (4 ** len(cut), 4 ** len(other))
    realigned = matrix.reshape([2] * 8).transpose(axes).reshape(shape)
    singular = np.linalg.svd(realigned, compute_uv=False)
    return int(np.count_nonzero(singular > 1e-9)), float(
        singular[singular > 1e-9][-1])


def feasible_graphs(bounds):
    return [graph for graph in GRAPHS
            if all(sum((a in cut) != (b in cut) for a, b in graph) >= bound
                   for cut, bound in zip(CUTS, bounds))]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    split = 13 + len(exact_block()) + len(reverse_block())
    prefix, suffix = source["gates"][:split], source["gates"][split:]
    assert sum(gate["gate"] == "cx" for gate in prefix) == 10
    assert sum(gate["gate"] == "cx" for gate in suffix) == 5
    s = unitary(suffix)
    records = []
    for a, b in itertools.combinations(range(16), 2):
        if labels[a] != labels[b] or (a ^ b).bit_count() < 2:
            continue
        for angle_index in (1, 2):
            angle = angle_index * math.pi / 4
            for phase_index in (0, 1):
                phase = phase_index * math.pi / 2
                gauged = two_level(a, b, angle, phase) @ s
                ranks = [schmidt_rank(gauged, cut) for cut in CUTS]
                bounds = [int(math.ceil(math.log2(rank)))
                          for rank, _ in ranks]
                graphs = feasible_graphs(bounds)
                records.append({"row_pair": [a, b], "eigenvalue_label": labels[a],
                                "hamming_distance": (a ^ b).bit_count(),
                                "angle_over_pi": angle_index / 4,
                                "phase_over_pi": phase_index / 2,
                                "ranks": [rank for rank, _ in ranks],
                                "smallest_nonzero_singular": [value for _, value in ranks],
                                "cut_cnot_bounds": bounds,
                                "feasible_graph_count": len(graphs),
                                "example_graphs": graphs[:8]})
    baseline_ranks = [schmidt_rank(s, cut)[0] for cut in CUTS]
    result = {"source": "topology16_15_exact_candidate.json",
              "split_gate": split, "prefix_cnot": 10, "suffix_cnot": 5,
              "cut_order": [list(cut) for cut in CUTS],
              "baseline_suffix_ranks": baseline_ranks,
              "gauge_trials": len(records),
              "rank_compatible_trials": sum(r["feasible_graph_count"] > 0
                                            for r in records),
              "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))
    print(json.dumps({"top_compatible": sorted(records,
                                                key=lambda r: -r["feasible_graph_count"])[:5]},
                     indent=2))


if __name__ == "__main__":
    main()
