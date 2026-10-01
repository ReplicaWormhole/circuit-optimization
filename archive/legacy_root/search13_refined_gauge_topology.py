"""Numerical cut-topology screen for the refined simultaneous-rank gauge.

The run297 gauge is numerical, so its 127 ranks at threshold 1e-9 are
candidate-selection data, not exact certificates. We nonetheless apply
the exact CNOT-factor graph mincut to find seven-CNOT pair schedules that
could synthesize that candidate tail after the fixed six-CNOT prefix.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

from search13_fulltail_rank8_topology import BASE, reverse_cone
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_mixedcut_gauge_refine_result.json"
OUT = ROOT / "search13_refined_gauge_topology_result.json"
EDGES = tuple(combinations(range(4), 2))


def distance_to_one_deletion(schedule):
    return min(sum(a != b for a, b in zip(schedule, BASE[:k] + BASE[k+1:]))
               for k in range(8))


def main():
    data = json.loads(SOURCE.read_text())
    ranks = {int(k): int(v) for k, v in
             data["numerical_rank_profiles_at_1e_minus_9"].items()}
    assert len(ranks) == 127 and sorted(set(ranks.values())) == [2, 4, 8, 16]
    masks = sorted(ranks, key=lambda mask: (-ranks[mask], mask))
    forbidden = {tuple(row) for row in json.loads((ROOT / "search13_fulltail_invariant_commutant_result.json").read_text())
                 ["certified_no_commutant_subsets"]}
    counts = Counter()
    first_fail = Counter()
    survivors = []
    distances = Counter()
    for schedule in product(EDGES, repeat=7):
        if any(reverse_cone(schedule, j) in forbidden for j in range(4)):
            counts["exact_commutant_rejected"] += 1
            continue
        graph = bond_graph(schedule)
        for mask in masks:
            if ranks[mask] > 2**edge_connectivity_up_to_four(graph, mask):
                counts["numerical_cut_rejected"] += 1
                first_fail[mask] += 1
                break
        else:
            counts["survives_numerical_cut_screen"] += 1
            survivors.append([list(edge) for edge in schedule])
            distances[distance_to_one_deletion(schedule)] += 1
    assert sum(counts.values()) == 6**7
    result = {"source": SOURCE.name,
              "all_ordered_unoriented_seven_pair_schedules": 6**7,
              "numerical_rank_threshold": 1e-9,
              "target_rank_histogram": dict(sorted(Counter(ranks.values()).items())),
              "exclusive_counts": dict(counts),
              "first_failing_mask_counts": {str(k): v for k, v in sorted(first_fail.items())},
              "survivor_distance_to_one_deletion": dict(sorted(distances.items())),
              "surviving_schedules": survivors,
              "scope": "numerical gauge-rank candidate selection; cut rejection not rigorous until gauge exactified"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"counts": result["exclusive_counts"],
                      "distances": result["survivor_distance_to_one_deletion"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
