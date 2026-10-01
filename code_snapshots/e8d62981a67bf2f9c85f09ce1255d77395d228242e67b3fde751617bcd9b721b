"""Exact necessary seven-CNOT topology screen for the rank-eight gauged tail.

Operator-Schmidt ranks 4,4,4,4,14,14,8 impose crossing capacities;
the fixed-prefix commutant screen imposes reverse support-cone conditions.
The finite enumeration is orientation-independent and does not synthesize.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGET = ROOT / "search13_fulltail_invariant_exact_witness_result.json"
COMMUTANT = ROOT / "search13_fulltail_invariant_commutant_result.json"
OUT = ROOT / "search13_fulltail_rank8_topology_result.json"
EDGES = tuple(combinations(range(4), 2))
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
BASE = ((0, 2), (0, 2), (1, 3), (1, 3),
        (0, 1), (0, 1), (2, 3), (2, 3))


def crossing(schedule, cut):
    half = set(cut)
    return sum((a in half) != (b in half) for a, b in schedule)


def reverse_cone(schedule, wire):
    support = {wire}
    for edge in reversed(schedule):
        if support.intersection(edge):
            support.update(edge)
    return tuple(sorted(support))


def edit_distance_to_one_deletion(schedule):
    return min(sum(a != b for a, b in zip(schedule, BASE[:k] + BASE[k+1:]))
               for k in range(8))


def main():
    target = json.loads(TARGET.read_text())
    required = tuple(target["operator_schmidt_ranks"]["".join(map(str, cut))]
                     for cut in CUTS)
    forbidden = {tuple(s) for s in json.loads(COMMUTANT.read_text())
                 ["certified_no_commutant_subsets"]}
    stage = Counter()
    profiles = Counter()
    best = {}
    for schedule in product(EDGES, repeat=7):
        crosses = tuple(crossing(schedule, cut) for cut in CUTS)
        if any(rank > 2**count for rank, count in zip(required, crosses)):
            stage["cut_rank_rejected"] += 1
            continue
        if any(reverse_cone(schedule, j) in forbidden for j in range(4)):
            stage["commutant_cone_rejected"] += 1
            continue
        stage["survives_both"] += 1
        profile = crosses[4:]
        profiles["-".join(map(str, profile))] += 1
        degree_spread = max(crosses[:4]) - min(crosses[:4])
        repeated_adjacent = sum(schedule[i] == schedule[i+1] for i in range(6))
        score = (edit_distance_to_one_deletion(schedule),
                 degree_spread, repeated_adjacent, schedule)
        if profile not in best or score < best[profile][0]:
            best[profile] = (score, schedule, crosses)
    assert sum(stage.values()) == 6**7
    representatives = [
        {"balanced_crossings": list(profile),
         "single_wire_degrees": list(item[2][:4]),
         "edit_distance_to_exact8_one_deletion": item[0][0],
         "schedule": [list(edge) for edge in item[1]]}
        for profile, item in sorted(best.items())]
    result = {
        "target_operator_schmidt_ranks": required,
        "cuts": [list(cut) for cut in CUTS],
        "all_seven_cnot_pair_schedules": 6**7,
        "exclusive_counts": dict(stage),
        "survivor_balanced_crossing_profiles": dict(sorted(profiles.items())),
        "closest_representative_by_profile": representatives,
        "scope": "exact necessary capacity and support-cone screen; no gate synthesis",
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"counts": dict(stage), "profiles": dict(profiles)},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
