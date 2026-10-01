"""All 5460 disjoint two-row-swap output gauges near one exact-tail deletion.

Each output-row permutation preserves cycle diagonalization. Modular mixed
cut ranks give rigorous necessary conditions for the 124 seven-CNOT
one-deletion/one-rewire pair schedules, independent of local gate angles.
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
OUT = ROOT / "search13_double_row_swap_near_tail_rank_result.json"
PRIMES = (97, 193)


def disjoint_pairings():
    for a, b, c, d in combinations(range(16), 4):
        yield ((a, b), (c, d))
        yield ((a, c), (b, d))
        yield ((a, d), (b, c))


def main():
    assert not OUT.exists(), OUT
    bases = {p: target_mod(p) for p in PRIMES}
    near = schedules()
    assert len(near) == 124
    graphs = [bond_graph(s) for s in near]
    caps = {m: np.array([edge_connectivity_up_to_four(g, m) for g in graphs],
                        dtype=np.uint8) for m in MASKS}
    masks = sorted(MASKS, key=lambda m: (sum(caps[m] >= 4), int(caps[m].sum()), m))
    rows = []
    checks = Counter()
    count = 0
    for pairs in disjoint_pairings():
        count += 1
        permutation = list(range(16))
        for a, b in pairs:
            permutation[a], permutation[b] = permutation[b], permutation[a]
        possible = np.ones(len(near), dtype=bool)
        for prime in PRIMES:
            if not possible.any():
                break
            target = bases[prime][permutation, :]
            for mask in masks:
                rank = rank_mod(realignment(target, mask), prime)
                checks[prime] += 1
                possible &= rank <= (1 << caps[mask])
                if not possible.any():
                    break
        if possible.any():
            rows.append({"swapped_output_row_pairs": pairs,
                         "surviving_schedules": [[list(edge) for edge in near[i]]
                                                 for i in np.flatnonzero(possible)]})
        if count % 1000 == 0:
            print(json.dumps({"gauges_checked": count,
                              "gauges_with_survivors": len(rows)}), flush=True)
    assert count == 5460
    result = {"source": "exact14 eight-CNOT tail under two disjoint output-row swaps",
              "primes": PRIMES, "masks": len(MASKS),
              "row_swap_gauges": count, "neighborhood_schedules": len(near),
              "rank_checks_by_prime": dict(checks),
              "gauges_with_survivors": len(rows),
              "total_survivor_pairs": sum(len(r["surviving_schedules"]) for r in rows),
              "rows_with_survivors": rows,
              "scope": "restricted exact modular necessary-condition screen; no synthesis or global no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gauges_with_survivors": len(rows),
                      "total_survivor_pairs": result["total_survivor_pairs"],
                      "rank_checks": dict(checks)}), flush=True)


if __name__ == "__main__":
    main()
