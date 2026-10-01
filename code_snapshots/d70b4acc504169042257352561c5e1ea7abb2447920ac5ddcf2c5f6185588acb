"""Exact modular near-tail screen for three-output-CNOT eigenrow gauges.

Deduplicates the noncancelling directed three-CNOT sequences by their
16-row permutation before testing 124 near-tail seven-CNOT schedules.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

from search13_output_cx_near_tail_rank import EDGES, MASKS, PRIMES, output_cx_rows, target_mod
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment, schedules
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_three_output_cx_near_tail_rank_result.json"


def unique_gauges():
    gauges = {}
    sequence_count = 0
    for sequence in product(EDGES, repeat=3):
        if sequence[0] == sequence[1] or sequence[1] == sequence[2]:
            continue
        sequence_count += 1
        rows = output_cx_rows(*sequence[0])
        for edge in sequence[1:]:
            rows = rows[output_cx_rows(*edge)]
        gauges.setdefault(tuple(map(int, rows)), sequence)
    assert sequence_count == 1452
    return gauges, sequence_count


def main():
    assert not OUT.exists(), OUT
    gauges, sequence_count = unique_gauges()
    bases = {p: target_mod(p) for p in PRIMES}
    nearby = schedules()
    assert len(nearby) == 124
    graphs = {s: bond_graph(s) for s in nearby}
    capacities = {s: {m: edge_connectivity_up_to_four(graph, m)
                      for m in MASKS} for s, graph in graphs.items()}
    rows = []
    for number, (permutation, sequence) in enumerate(gauges.items(), 1):
        ranks = {}
        for p in PRIMES:
            target = bases[p][list(permutation), :]
            ranks[p] = {m: rank_mod(realignment(target, m), p) for m in MASKS}
        lower = {m: max(ranks[p][m] for p in PRIMES) for m in MASKS}
        survivors = [[list(edge) for edge in s] for s in nearby
                     if all(lower[m] <= 2**capacities[s][m] for m in MASKS)]
        rows.append({"representative_sequence": sequence,
                     "row_permutation": permutation,
                     "surviving_schedules": survivors,
                     "mod97_rank_histogram": dict(Counter(ranks[97].values())),
                     "mod193_rank_histogram": dict(Counter(ranks[193].values()))})
        if number % 100 == 0:
            print(json.dumps({"unique_gauges_checked": number,
                              "gauges_with_survivors": sum(bool(r["surviving_schedules"]) for r in rows)}), flush=True)
    result = {"source": "exact14 eight-CNOT tail followed by three output CNOTs",
              "primes": PRIMES, "masks": len(MASKS),
              "noncancelling_sequences": sequence_count,
              "unique_row_permutations": len(gauges),
              "neighborhood_schedules": len(nearby),
              "gauges_with_survivors": sum(bool(r["surviving_schedules"]) for r in rows),
              "total_survivor_pairs": sum(len(r["surviving_schedules"]) for r in rows),
              "rows": rows,
              "scope": "restricted exact modular necessary-condition screen; no general 13-CNOT no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"unique_gauges": len(gauges),
                      "gauges_with_survivors": result["gauges_with_survivors"],
                      "total_survivor_pairs": result["total_survivor_pairs"]}), flush=True)


if __name__ == "__main__":
    main()
