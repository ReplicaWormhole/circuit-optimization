"""Exact modular cut screen for one output CNOT followed by seven-tail reuse.

An output CNOT permutes eigenrows of the exact 14-CNOT diagonalizer. For
each of 12 directed output CNOTs, test the 124 seven-pair schedules made
by deleting one CNOT from its eight-CNOT tail and changing at most one
remaining pair. Modular ranks lower-bound exact ranks; CNOT-factor mincuts
upper-bound any circuit with a fixed schedule and arbitrary local gates.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_rank8_neighborhood_mixedcut import (
    MASKS, PRIMES, field_mod, rank_mod, realignment, schedules,
)
from search13_rank8_global_twocut import bond_graph, edge_connectivity_up_to_four
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_output_cx_near_tail_rank_result.json"
EDGES = tuple((c, t) for c in range(4) for t in range(4) if c != t)


def output_cx_rows(control, target):
    rows = np.arange(16, dtype=np.int64)
    flip = 1 << (3 - target)
    active = (rows >> (3 - control)) & 1
    return rows ^ (active * flip)


def target_mod(prime):
    field = CyclotomicPair()
    tail = tail_matrix(field)
    root = pow(primitive_root(prime), (prime - 1) // 96, prime)
    assert pow(root, 96, prime) == 1 and pow(root, 48, prime) != 1
    return np.array([[field_mod(v, root, prime) for v in row]
                     for row in tail], dtype=np.int64)


def profile_for_output_cx(base, edge, prime):
    matrix = base[output_cx_rows(*edge), :]
    return {mask: rank_mod(realignment(matrix, mask), prime)
            for mask in MASKS}


def main():
    assert not OUT.exists(), OUT
    bases = {p: target_mod(p) for p in PRIMES}
    nearby = schedules()
    assert len(nearby) == 124
    graphs = {schedule: bond_graph(schedule) for schedule in nearby}
    mincuts = {schedule: {mask: edge_connectivity_up_to_four(graph, mask)
                          for mask in MASKS}
               for schedule, graph in graphs.items()}
    rows = []
    for edge in EDGES:
        by_prime = {p: profile_for_output_cx(bases[p], edge, p) for p in PRIMES}
        lower = {mask: max(by_prime[p][mask] for p in PRIMES) for mask in MASKS}
        counts = Counter()
        survivors = []
        for schedule in nearby:
            first = next((mask for mask in MASKS
                          if lower[mask] > 2**mincuts[schedule][mask]), None)
            if first is None:
                survivors.append([list(pair) for pair in schedule])
                counts["passes_modular_cut_tests"] += 1
            else:
                counts["modular_cut_rejected"] += 1
        assert sum(counts.values()) == len(nearby)
        row = {"output_cx": edge, "counts": dict(counts),
               "rank_histogram_mod97": dict(Counter(by_prime[97].values())),
               "rank_histogram_mod193": dict(Counter(by_prime[193].values())),
               "surviving_schedules": survivors}
        rows.append(row)
        print(json.dumps({"output_cx": edge, "counts": counts}), flush=True)
    result = {"source": "exact 14-CNOT tail with one output CNOT row permutation",
              "primes": PRIMES, "masks": len(MASKS),
              "neighborhood_schedules": len(nearby), "rows": rows,
              "scope": "exact necessary modular-rank screen for one output CNOT and 124 near-tail pair orders; survivors require further validation"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"total_survivors": sum(len(row["surviving_schedules"]) for row in rows),
                      "gauges_with_survivors": sum(bool(row["surviving_schedules"]) for row in rows)}), flush=True)


if __name__ == "__main__":
    main()
