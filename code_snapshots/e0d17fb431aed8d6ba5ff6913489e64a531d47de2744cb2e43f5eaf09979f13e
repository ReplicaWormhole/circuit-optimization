"""All output-row three-cycles versus seven-tail orders within two rewires.

This extends the exact necessary screen of run 322 from 124 near-tail orders
to the 1632 orders one deletion and at most two pair rewires from the exact
eight-CNOT tail. A surviving gauge/order is not a synthesized circuit.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np

from search13_all_linear_output_distance2 import neighborhood
from search13_output_cx_near_tail_rank import MASKS, target_mod
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment
from search13_threecycle_near_tail_rank import three_cycles


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_threecycle_distance2_rank_result.json"
PRIMES = (97, 193)


def main():
    assert not OUT.exists(), OUT
    near = neighborhood()
    assert len(near) == 1632
    graphs = [bond_graph(s) for s, _ in near]
    caps = {m: np.array([edge_connectivity_up_to_four(g, m) for g in graphs],
                        dtype=np.uint8) for m in MASKS}
    masks = sorted(MASKS, key=lambda m: (sum(caps[m] >= 4),
                                       int(caps[m].sum()), m))
    bases = {p: target_mod(p) for p in PRIMES}
    counts = Counter()
    checks = Counter()
    survivors = []
    for cycle, permutation in three_cycles():
        possible = np.ones(len(near), dtype=bool)
        for prime in PRIMES:
            if not possible.any():
                break
            target = bases[prime][list(permutation), :]
            for mask in masks:
                rank = rank_mod(realignment(target, mask), prime)
                checks[prime] += 1
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
    distinct = sorted({i for row in survivors
                       for i in row["surviving_schedule_indices"]})
    result = {"source": "exact 14-CNOT eight-CNOT tail",
              "primes": PRIMES, "masks": len(MASKS),
              "three_cycle_gauges": 1120,
              "neighborhood_schedules": len(near),
              "schedule_order": [[list(edge) for edge in s] for s, _ in near],
              "schedule_distances": [d for _, d in near],
              "rank_checks_by_prime": dict(checks),
              "survivors_by_gauge": dict(sorted(counts.items())),
              "gauges_with_survivors": len(survivors),
              "total_survivor_pairs": sum(len(row["surviving_schedule_indices"])
                                          for row in survivors),
              "distinct_surviving_schedule_indices": distinct,
              "rows_with_survivors": survivors,
              "scope": "exact necessary modular mixed-cut screen of finite row-gauge/topology family; no synthesis or global 13-CNOT bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gauges_with_survivors": len(survivors),
                      "total_survivor_pairs": result["total_survivor_pairs"],
                      "distinct_surviving_schedules": len(distinct),
                      "rank_checks": dict(checks)}), flush=True)


if __name__ == "__main__":
    main()
