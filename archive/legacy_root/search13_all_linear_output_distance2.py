"""All linear output gauges against seven-tail orders within two rewires.

The set comprises seven-CNOT pair orders obtainable by deleting one gate
from the exact eight-CNOT tail and changing at most two retained pairs.
Screen the 20160 GL(4,2) output row permutations by modular mixed-cut
ranks and CNOT-factor mincuts. Survival is only a necessary condition.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from search13_all_linear_output_near_tail_rank import linear_row_permutations
from search13_fulltail_rank8_topology import BASE
from search13_output_cx_near_tail_rank import MASKS, target_mod
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_all_linear_output_distance2_result.json"
EDGES = tuple(combinations(range(4), 2))
PRIMES = (97, 193)


def neighborhood():
    deleted = tuple(BASE[:k] + BASE[k + 1:] for k in range(8))
    schedules = []
    for s in product(EDGES, repeat=7):
        distance = min(sum(a != b for a, b in zip(s, old)) for old in deleted)
        if distance <= 2:
            schedules.append((s, distance))
    assert len(schedules) == 1632
    assert sum(d <= 1 for _, d in schedules) == 124
    return schedules


def main():
    assert not OUT.exists(), OUT
    gauges = linear_row_permutations()
    near = neighborhood()
    graphs = [bond_graph(schedule) for schedule, _ in near]
    caps = {m: np.array([edge_connectivity_up_to_four(g, m) for g in graphs],
                        dtype=np.uint8) for m in MASKS}
    masks = sorted(MASKS, key=lambda m: (sum(caps[m] >= 4),
                                       int(caps[m].sum()), m))
    base = {p: target_mod(p) for p in PRIMES}
    hist = Counter()
    checks = Counter()
    rows = []
    for index, (permutation, (depth, word)) in enumerate(gauges.items(), 1):
        possible = np.ones(len(near), dtype=bool)
        for prime in PRIMES:
            if not possible.any():
                break
            target = base[prime][list(permutation), :]
            for mask in masks:
                rank = rank_mod(realignment(target, mask), prime)
                checks[prime] += 1
                possible &= rank <= (1 << caps[mask])
                if not possible.any():
                    break
        indices = np.flatnonzero(possible)
        hist[(depth, len(indices))] += 1
        if len(indices):
            rows.append({"row_permutation": permutation,
                         "minimum_output_cx_length": depth,
                         "representative_output_cxs": word,
                         "surviving_indices": indices.tolist(),
                         "distance_histogram": dict(Counter(near[i][1] for i in indices))})
        if index % 2000 == 0:
            print(json.dumps({"gauges_checked": index,
                              "gauges_with_survivors": len(rows)}), flush=True)
    result = {"source": "exact14 eight-CNOT tail, GL(4,2) output row gauges",
              "primes": PRIMES, "masks": len(MASKS),
              "linear_output_gauges": len(gauges),
              "neighborhood": "one deletion plus at most two retained-pair rewires",
              "neighborhood_schedules": len(near),
              "schedule_order": [[list(edge) for edge in s] for s, _ in near],
              "schedule_distances": [d for _, d in near],
              "rank_checks_by_prime": dict(checks),
              "survival_histogram_depth_and_count": {f"{d},{s}": c for (d, s), c in sorted(hist.items())},
              "gauges_with_survivors": len(rows),
              "total_survivor_pairs": sum(len(r["surviving_indices"]) for r in rows),
              "rows_with_survivors": rows,
              "scope": "exact necessary mixed-cut screen of finite output-gauge/topology family; no synthesis or global 13-CNOT no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gauges_with_survivors": len(rows),
                      "total_survivor_pairs": result["total_survivor_pairs"],
                      "rank_checks": dict(checks)}), flush=True)


if __name__ == "__main__":
    main()
