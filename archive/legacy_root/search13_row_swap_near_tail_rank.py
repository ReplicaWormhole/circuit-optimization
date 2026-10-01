"""All 120 exact output-row transpositions against near-tail seven-CNOTs.

Permuting output rows preserves diagonalization (with permuted eigenvalue
labels). Test every single row transposition of the exact14 eight-CNOT tail
against the 124 one-deletion/one-rewire seven-CNOT pair schedules. Modular
rank lower bounds and CNOT-factor mincuts are necessary conditions only.
"""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

from search13_output_cx_near_tail_rank import MASKS, PRIMES, target_mod
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment, schedules
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_row_swap_near_tail_rank_result.json"


def main():
    assert not OUT.exists(), OUT
    bases = {p: target_mod(p) for p in PRIMES}
    near = schedules()
    assert len(near) == 124
    graphs = {s: bond_graph(s) for s in near}
    caps = {s: {m: edge_connectivity_up_to_four(g, m) for m in MASKS}
            for s, g in graphs.items()}
    rows = []
    for a, b in combinations(range(16), 2):
        permutation = list(range(16))
        permutation[a], permutation[b] = permutation[b], permutation[a]
        rank = {}
        for p in PRIMES:
            target = bases[p][permutation, :]
            rank[p] = {m: rank_mod(realignment(target, m), p) for m in MASKS}
        lower = {m: max(rank[p][m] for p in PRIMES) for m in MASKS}
        survivors = [[list(edge) for edge in s] for s in near
                     if all(lower[m] <= 2**caps[s][m] for m in MASKS)]
        rows.append({"swapped_output_rows": (a, b),
                     "surviving_schedules": survivors,
                     "rank_histogram_mod97": dict(Counter(rank[97].values())),
                     "rank_histogram_mod193": dict(Counter(rank[193].values()))})
    result = {"source": "exact14 eight-CNOT tail under all 120 output-row transpositions",
              "primes": PRIMES, "masks": len(MASKS),
              "row_swap_gauges": len(rows), "neighborhood_schedules": len(near),
              "gauges_with_survivors": sum(bool(r["surviving_schedules"]) for r in rows),
              "total_survivor_pairs": sum(len(r["surviving_schedules"]) for r in rows),
              "rows": rows,
              "scope": "restricted exact modular necessary-condition screen; no synthesis or global no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gauges_with_survivors": result["gauges_with_survivors"],
                      "total_survivor_pairs": result["total_survivor_pairs"]}), flush=True)


if __name__ == "__main__":
    main()
