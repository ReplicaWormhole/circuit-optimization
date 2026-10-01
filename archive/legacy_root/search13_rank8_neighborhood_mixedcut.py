"""Exact modular mixed-cut screen near one-deletion eight-CNOT tails.

The target is the exact rank-eight gauged tail. We enumerate every ordered
seven-CNOT pair schedule obtained by deleting one gate from its original
eight-CNOT tail and changing at most one remaining pair. For every mixed
input/output bipartition, a finite-field target realignment rank is a
rigorous lower bound on its characteristic-zero rank. A tensor-network
wire-bond mincut gives a rigorous upper bound for an arbitrary-local-gate
circuit with that ordered CNOT pair schedule and either directions.
"""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_fulltail_invariant_exact_witness import gauge_matrix
from search13_fulltail_rank8_topology import (
    BASE, CUTS, crossing, reverse_cone,
)
from search13_fulltail_rank8_mixedcut import wire_graph
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_rank8_neighborhood_mixedcut_result.json"
EDGES = tuple(combinations(range(4), 2))
PRIMES = (97, 193)
MASKS = tuple(range(1, 255, 2))


def field_mod(value, root, prime):
    out = 0
    for coefficient in value.rep:
        numerator = int(coefficient.numerator) % prime
        denominator = int(coefficient.denominator) % prime
        assert denominator
        out = (out * root + numerator * pow(denominator, -1, prime)) % prime
    return out


def rank_mod(matrix, prime):
    a = np.array(matrix, dtype=np.int64, copy=True) % prime
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        candidates = np.flatnonzero(a[rank:, col])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, prime) % prime
        if rank + 1 < rows:
            factors = a[rank + 1:, col].copy()
            a[rank + 1:] = (a[rank + 1:] - factors[:, None] * a[rank]) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def realignment(matrix, mask):
    left = tuple(j for j in range(8) if mask & (1 << j))
    right = tuple(j for j in range(8) if j not in left)
    return matrix.reshape((2,) * 8).transpose(left + right).reshape(
        (1 << len(left), 1 << len(right)))


def target_rank_lower_bounds():
    exact = CyclotomicPair()
    gauge, _ = gauge_matrix(exact)
    target = exact.matmul(gauge, tail_matrix(exact))
    ranks = {mask: 0 for mask in MASKS}
    for prime in PRIMES:
        root = pow(primitive_root(prime), (prime - 1) // 96, prime)
        assert pow(root, 96, prime) == 1 and pow(root, 48, prime) != 1
        numeric = np.array([[field_mod(v, root, prime) for v in row]
                            for row in target], dtype=np.int64)
        for mask in MASKS:
            ranks[mask] = max(ranks[mask], rank_mod(realignment(numeric, mask), prime))
    return ranks


def schedules():
    result = set()
    for omitted in range(8):
        base = BASE[:omitted] + BASE[omitted + 1:]
        result.add(base)
        for index in range(7):
            for edge in EDGES:
                if edge == base[index]:
                    continue
                altered = list(base)
                altered[index] = edge
                result.add(tuple(altered))
    return tuple(sorted(result))


def mincuts(edges):
    boundary = ((np.array(MASKS, dtype=np.uint16)[:, None]
                 >> np.arange(8, dtype=np.uint16)) & 1).astype(np.int8)
    gate = ((np.arange(128, dtype=np.uint16)[:, None]
             >> np.arange(7, dtype=np.uint16)) & 1).astype(np.int8)
    sides = np.concatenate((np.broadcast_to(boundary[:, None, :], (127, 128, 8)),
                            np.broadcast_to(gate[None, :, :], (127, 128, 7))),
                           axis=2)
    cuts = np.zeros((127, 128), dtype=np.int8)
    for a, b in edges:
        cuts += sides[:, :, a] != sides[:, :, b]
    return cuts.min(axis=1)


def main():
    ranks = target_rank_lower_bounds()
    target = json.loads((ROOT / "search13_fulltail_invariant_exact_witness_result.json").read_text())
    required = tuple(target["operator_schmidt_ranks"]["".join(map(str, cut))]
                     for cut in CUTS)
    commutant = json.loads((ROOT / "search13_fulltail_invariant_commutant_result.json").read_text())
    forbidden = {tuple(row) for row in commutant["certified_no_commutant_subsets"]}
    rows = []
    survivors = []
    stage = {"neighborhood": 0, "ordinary_cut_rejected": 0,
             "commutant_rejected": 0, "mixed_cut_rejected": 0,
             "survives_all": 0}
    for schedule in schedules():
        stage["neighborhood"] += 1
        crosses = tuple(crossing(schedule, cut) for cut in CUTS)
        if any(rank > 2**count for rank, count in zip(required, crosses)):
            stage["ordinary_cut_rejected"] += 1
            continue
        if any(reverse_cone(schedule, j) in forbidden for j in range(4)):
            stage["commutant_rejected"] += 1
            continue
        caps = mincuts(wire_graph(schedule))
        violations = [(mask, ranks[mask], int(cap))
                      for mask, cap in zip(MASKS, caps)
                      if ranks[mask] > 2**int(cap)]
        if violations:
            stage["mixed_cut_rejected"] += 1
            rows.append({"schedule": [list(edge) for edge in schedule],
                         "first_violating_cut": list(violations[0]),
                         "violating_cut_count": len(violations)})
        else:
            stage["survives_all"] += 1
            survivors.append([list(edge) for edge in schedule])
    assert stage["neighborhood"] == sum(v for k, v in stage.items() if k != "neighborhood")
    result = {"source": "exact rank-eight output gauge and exact eight-CNOT tail",
              "primes": PRIMES, "source_eight_pair_schedule": BASE,
              "neighborhood": "one deletion plus at most one changed pair",
              "counts": stage,
              "target_modular_rank_lower_bounds": {str(k): v for k, v in ranks.items()},
              "mixed_cut_exclusions": rows,
              "surviving_schedules": survivors,
              "scope": "exact necessary cut screen for listed neighborhood schedules only; modular rank is a lower bound, mincut an upper bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stage, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
