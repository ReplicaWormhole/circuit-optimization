"""Exhaust all within-eigenspace row permutations for one mixed-cut rank.

The exact target lies in Q(zeta_96). Reduction modulo 97 provides a lower
bound on its characteristic-zero rank. A rank greater than eight at even one
mixed cut rules out a rank-eight target for that row-permutation gauge.
"""

from collections import Counter
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
from numba import njit

from search13_fulltail_invariant_gauge_rank import setup
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment
from search13_sparse_gauge_clifford_pairs import target_mod


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_sector_row_permutations_exhaustive_result.json"
P = 97
MASK = 83


@njit(cache=True)
def rank_realign_permuted(base, perm, indices, prime):
    a = np.empty((16, 16), dtype=np.int64)
    for r in range(16):
        for c in range(16):
            k = indices[r, c]
            a[r, c] = base[perm[k // 16], k % 16]
    rank = 0
    for col in range(16):
        pivot = -1
        for row in range(rank, 16):
            if a[row, col] != 0:
                pivot = row
                break
        if pivot < 0:
            continue
        if pivot != rank:
            for j in range(col, 16):
                a[rank, j], a[pivot, j] = a[pivot, j], a[rank, j]
        inv = pow_mod(a[rank, col], prime - 2, prime)
        for j in range(col, 16):
            a[rank, j] = (a[rank, j] * inv) % prime
        for row in range(rank + 1, 16):
            factor = a[row, col]
            if factor != 0:
                for j in range(col, 16):
                    a[row, j] = (a[row, j] - factor * a[rank, j]) % prime
        rank += 1
        if rank == 16:
            break
    return rank


@njit(cache=True)
def pow_mod(base, exponent, prime):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = answer * base % prime
        base = base * base % prime
        exponent //= 2
    return answer


def index_table(mask):
    left = tuple(j for j in range(8) if mask & (1 << j))
    right = tuple(j for j in range(8) if j not in left)
    return np.arange(256).reshape((2,) * 8).transpose(left + right).reshape(16, 16)


def main():
    assert not OUT.exists(), OUT
    base, _, _ = target_mod(P)
    _, _, _, sectors = setup()
    indices = index_table(MASK)
    per_sector = [list(permutations(s)) for s in sectors]
    assert np.prod([len(x) for x in per_sector]) == 622080
    perm = np.arange(16)
    # Check the compiled reducer against the independent Python reducer.
    for choices in (tuple(s for s in sectors), tuple(x[-1] for x in per_sector)):
        for sector, assigned in zip(sectors, choices):
            perm[np.asarray(sector)] = assigned
        assert rank_realign_permuted(base, perm, indices, P) == rank_mod(realignment(base[perm, :], MASK), P)
    counts = Counter()
    minimum = (99, None)
    for choices in product(*per_sector):
        for sector, assigned in zip(sectors, choices):
            perm[np.asarray(sector)] = assigned
        rank = int(rank_realign_permuted(base, perm, indices, P))
        counts[rank] += 1
        if rank < minimum[0]:
            minimum = (rank, perm.tolist())
            print(json.dumps({"new_minimum": rank, "seen": sum(counts.values())}), flush=True)
    result = {"source": "exact run266 gauge times exact14 tail",
              "family": "all independent output row permutations within shift eigenspaces",
              "sector_sizes": list(map(len, sectors)),
              "total": sum(counts.values()), "cut_mask": MASK, "prime": P,
              "rank_histogram_mod97": dict(sorted(counts.items())),
              "minimum_mod97_rank": minimum[0],
              "minimum_example_row_permutation": minimum[1],
              "scope": "exhaustive finite output permutation family; modular rank is a rigorous characteristic-zero lower bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"minimum": minimum[0], "total": result["total"]}), flush=True)


if __name__ == "__main__":
    main()
