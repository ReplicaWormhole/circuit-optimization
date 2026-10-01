"""Exhaust all seven-CNOT pair schedules against 127 exact mixed cuts.

For the fixed exact rank-eight-gauged tail, modular realignment ranks
at 97 and 193 give rigorous characteristic-zero lower bounds. The
orientation-independent factor graph of CNOTs gives dimension-two
bond mincut upper bounds, with arbitrary local gates allowed.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

from search13_fulltail_rank8_topology import (
    CUTS, crossing, reverse_cone, edit_distance_to_one_deletion,
)
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_rank8_global_allcuts_result.json"
EDGES = tuple(combinations(range(4), 2))


def main():
    witness = json.loads((ROOT / "search13_fulltail_invariant_exact_witness_result.json").read_text())
    ordinary = tuple(witness["operator_schmidt_ranks"]["".join(map(str, cut))]
                     for cut in CUTS)
    ranks = {int(k): v for k, v in json.loads((ROOT / "search13_rank8_neighborhood_mixedcut_result.json").read_text())
             ["target_modular_rank_lower_bounds"].items()}
    masks = sorted(ranks, key=lambda mask: (-((ranks[mask]-1).bit_length()),
                                           -ranks[mask], mask))
    forbidden = {tuple(row) for row in json.loads((ROOT / "search13_fulltail_invariant_commutant_result.json").read_text())
                 ["certified_no_commutant_subsets"]}
    stage = Counter()
    first_cut = Counter()
    survivors = []
    survivor_distances = Counter()
    for schedule in product(EDGES, repeat=7):
        crosses = tuple(crossing(schedule, cut) for cut in CUTS)
        if any(rank > 2**count for rank, count in zip(ordinary, crosses)):
            stage["ordinary_cut_rejected"] += 1
            continue
        if any(reverse_cone(schedule, j) in forbidden for j in range(4)):
            stage["commutant_rejected"] += 1
            continue
        graph = bond_graph(schedule)
        for mask in masks:
            bonds = edge_connectivity_up_to_four(graph, mask)
            if ranks[mask] > 2**bonds:
                stage["mixed_cut_rejected"] += 1
                first_cut[mask] += 1
                break
        else:
            stage["survives_all"] += 1
            survivors.append([list(edge) for edge in schedule])
            survivor_distances[edit_distance_to_one_deletion(schedule)] += 1
    assert sum(stage.values()) == 6**7
    result = {"target": "exact rank-eight gauged tail after fixed six-CNOT prefix",
              "all_ordered_unoriented_seven_pair_schedules": 6**7,
              "mixed_boundary_masks_tested": len(masks),
              "mask_test_order": masks,
              "exclusive_counts": dict(stage),
              "first_violating_mask_counts": {str(k): v for k, v in sorted(first_cut.items())},
              "survivor_distance_to_one_deletion": dict(sorted(survivor_distances.items())),
              "surviving_schedules": survivors,
              "scope": "necessary exact finite-field rank and tensor mincut screens for this fixed exact gauge; no global 13-CNOT theorem"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["exclusive_counts"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
