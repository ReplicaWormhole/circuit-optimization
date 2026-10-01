"""Exact necessary screen for all output-row 3-cycles and near-tail orders.

Each of the 1120 oriented three-cycles permutes rows of the exact 14-CNOT
tail. For each gauge, modular mixed-cut ranks are lower bounds on exact
realignment ranks; CNOT-factor mincuts upper-bound the ranks attainable by
each of the 124 one-deletion/one-rewire seven-CNOT pair orders. Survival is
only a necessary condition for a fixed target and pair order.
"""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from search13_output_cx_near_tail_rank import MASKS, target_mod
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment, schedules
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_threecycle_near_tail_rank_result.json"
PRIME = 97


def three_cycles():
    for a, b, c in combinations(range(16), 3):
        for x, y, z in ((a, b, c), (a, c, b)):
            perm = list(range(16))
            perm[x], perm[y], perm[z] = y, z, x
            yield (x, y, z), tuple(perm)


def main():
    assert not OUT.exists(), OUT
    base = target_mod(PRIME)
    near = schedules()
    assert len(near) == 124
    graphs = [bond_graph(s) for s in near]
    caps = {m: np.array([edge_connectivity_up_to_four(g, m) for g in graphs],
                        dtype=np.uint8) for m in MASKS}
    masks = sorted(MASKS, key=lambda m: (sum(caps[m] >= 4),
                                       int(caps[m].sum()), m))
    counts = Counter()
    checks = 0
    survivors = []
    for cycle, permutation in three_cycles():
        target = base[list(permutation), :]
        possible = np.ones(len(near), dtype=bool)
        for mask in masks:
            rank = rank_mod(realignment(target, mask), PRIME)
            checks += 1
            possible &= rank <= (1 << caps[mask])
            if not possible.any():
                break
        indices = np.flatnonzero(possible).tolist()
        counts[len(indices)] += 1
        if indices:
            survivors.append({"output_cycle": cycle,
                              "row_permutation": permutation,
                              "surviving_schedule_indices": indices})
    assert sum(counts.values()) == 1120
    result = {"source": "exact 14-CNOT eight-CNOT tail",
              "prime": PRIME, "masks": len(MASKS),
              "three_cycle_gauges": 1120,
              "neighborhood_schedules": len(near),
              "schedule_order": [[list(edge) for edge in s] for s in near],
              "rank_checks": checks,
              "survivors_by_gauge": dict(sorted(counts.items())),
              "gauges_with_survivors": len(survivors),
              "total_survivor_pairs": sum(len(row["surviving_schedule_indices"])
                                          for row in survivors),
              "rows_with_survivors": survivors,
              "scope": "exact necessary modular mixed-cut screen of finite row-gauge/topology family; no global 13-CNOT bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gauges_with_survivors": len(survivors),
                      "total_survivor_pairs": result["total_survivor_pairs"],
                      "rank_checks": checks}), flush=True)


if __name__ == "__main__":
    main()
