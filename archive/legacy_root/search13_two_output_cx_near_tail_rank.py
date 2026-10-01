"""Modular near-tail screen for two noncancelling output CNOT gauges.

Applies each ordered pair of distinct directed output CNOTs to the exact
eight-CNOT tail, then tests the 124 one-deletion/one-rewire seven-CNOT pair
orders using exact modular-rank lower bounds and graph mincut upper bounds.
This is a restricted necessary-condition screen, not a global no-go.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

from search13_output_cx_near_tail_rank import (
    EDGES, PRIMES, MASKS, output_cx_rows, target_mod,
)
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment, schedules
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_two_output_cx_near_tail_rank_result.json"


def output_rows(sequence):
    first, second = sequence
    rows_first = output_cx_rows(*first)
    rows_second = output_cx_rows(*second)
    # (CX_second CX_first T)[y,:] = T[CX_first(CX_second(y)),:].
    return rows_first[rows_second]


def main():
    assert not OUT.exists(), OUT
    bases = {p: target_mod(p) for p in PRIMES}
    nearby = schedules()
    assert len(nearby) == 124
    # Every graph/mask cut is independent of output eigenrow gauge.
    graphs = {s: bond_graph(s) for s in nearby}
    mincuts = {s: {m: edge_connectivity_up_to_four(graph, m)
                   for m in MASKS} for s, graph in graphs.items()}
    sequences = tuple(pair for pair in product(EDGES, repeat=2)
                      if pair[0] != pair[1])
    assert len(sequences) == 132
    rows = []
    total = 0
    for sequence in sequences:
        permutation = output_rows(sequence)
        ranks = {}
        for p in PRIMES:
            matrix = bases[p][permutation, :]
            ranks[p] = {m: rank_mod(realignment(matrix, m), p) for m in MASKS}
        lower = {m: max(ranks[p][m] for p in PRIMES) for m in MASKS}
        survivors = []
        for s in nearby:
            if all(lower[m] <= 2**mincuts[s][m] for m in MASKS):
                survivors.append([list(edge) for edge in s])
        total += len(survivors)
        rows.append({"output_cxs": sequence,
                     "surviving_schedules": survivors,
                     "mod97_rank_histogram": dict(Counter(ranks[97].values())),
                     "mod193_rank_histogram": dict(Counter(ranks[193].values()))})
        if len(rows) % 12 == 0:
            print(json.dumps({"completed": len(rows), "survivors_so_far": total}), flush=True)
    result = {"source": "exact14 eight-CNOT tail followed by two output CNOTs",
              "primes": PRIMES, "masks": len(MASKS), "sequences": len(sequences),
              "neighborhood_schedules": len(nearby), "total_survivors": total,
              "gauges_with_survivors": sum(bool(r["surviving_schedules"]) for r in rows),
              "rows": rows,
              "scope": "restricted exact modular necessary-condition screen; survivor is not synthesis"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"total_survivors": total,
                      "gauges_with_survivors": result["gauges_with_survivors"]}), flush=True)


if __name__ == "__main__":
    main()
