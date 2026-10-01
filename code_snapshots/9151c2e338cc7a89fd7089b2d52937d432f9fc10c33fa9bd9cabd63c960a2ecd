"""Exact finite-field screen of sparse output-gauge row permutations.

Starting from the exact rank-eight gauge, apply two or three disjoint swaps
within V4 output eigenspaces. Test whether any of the eight seven-CNOT tail
orders excluded for the old gauge by directional mixed cuts becomes compatible
with the two decisive masks 83 and 163. Modular ranks are certified lower
bounds on characteristic-zero rank, not upper bounds.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_fulltail_invariant_exact_ansatz import LABELS, tail_matrix
from search13_fulltail_invariant_exact_witness import gauge_matrix
from search13_rank8_directional_mincut import graph, mincut
from search13_rank8_neighborhood_mixedcut import field_mod, rank_mod, realignment
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_rank8_neighborhood_mixedcut_result.json"
OUT = ROOT / "search13_sparse_gauge_swaps_result.json"
MASKS = (83, 163)
PRIMES = (97, 193)


def swaps():
    pairs = [(a, b) for a, b in combinations(range(16), 2)
             if LABELS[a] == LABELS[b]]
    families = []
    for size in (2, 3):
        families.extend(tuple(group) for group in combinations(pairs, size)
                        if len({row for pair in group for row in pair}) == 2 * size)
    assert len(families) == 1620
    return families


def permutation(group):
    perm = list(range(16))
    for a, b in group:
        perm[a], perm[b] = perm[b], perm[a]
    assert all(LABELS[row] == LABELS[old] for row, old in enumerate(perm))
    return perm


def targets_mod():
    exact = CyclotomicPair()
    g, _ = gauge_matrix(exact)
    target = exact.matmul(g, tail_matrix(exact))
    matrices = []
    for prime in PRIMES:
        root = pow(primitive_root(prime), (prime - 1) // 96, prime)
        assert pow(root, 96, prime) == 1 and pow(root, 48, prime) != 1
        matrices.append(np.asarray([[field_mod(v, root, prime) for v in row]
                                    for row in target], dtype=np.int64))
    return matrices


def orientation_caps(schedules):
    result = []
    for schedule in schedules:
        rows = []
        for directions in product(range(2), repeat=7):
            g = graph(schedule, directions)
            rows.append((tuple(directions), tuple(mincut(g, mask) for mask in MASKS)))
        result.append(rows)
    return result


def main():
    assert not OUT.exists()
    schedules = json.loads(SOURCE.read_text())["surviving_schedules"]
    assert len(schedules) == 8
    caps = orientation_caps(schedules)
    targets = targets_mod()
    base = tuple(max(rank_mod(realignment(target, mask), prime)
                     for target, prime in zip(targets, PRIMES)) for mask in MASKS)
    assert base == (14, 14)
    assert all(not any(all(base[i] <= 2**cuts[i] for i in range(2))
                       for _, cuts in rows) for rows in caps)

    count_by_size = Counter()
    rank_profiles = Counter()
    possible = []
    for group in swaps():
        count_by_size[len(group)] += 1
        perm = permutation(group)
        ranks = tuple(max(rank_mod(realignment(target[perm], mask), prime)
                          for target, prime in zip(targets, PRIMES))
                      for mask in MASKS)
        rank_profiles[ranks] += 1
        revived = []
        for index, orientations in enumerate(caps):
            admissible = [(direction, cuts) for direction, cuts in orientations
                          if all(ranks[i] <= 2**cuts[i] for i in range(2))]
            if admissible:
                revived.append({"schedule_index": index,
                                "direction_count": len(admissible),
                                "example_directions": admissible[0][0],
                                "example_mincuts": admissible[0][1]})
        if revived:
            possible.append({"swaps": group, "modular_rank_lower_bounds": ranks,
                             "revived_on_two_masks": revived})
    result = {"source": SOURCE.name,
              "base_gauge": "search13_fulltail_invariant_exact_witness.py",
              "family": "two or three disjoint same-eigenvalue row transpositions composed on the left of the exact rank-eight gauge",
              "field": "Q(zeta_96) reduced modulo 97 and 193",
              "primes": PRIMES, "masks": MASKS,
              "old_gauge_rank_lower_bounds": base,
              "family_counts": dict(count_by_size),
              "modular_rank_profile_counts": {f"{a},{b}": n for (a, b), n in sorted(rank_profiles.items())},
              "two_mask_possible": possible,
              "scope": "two-cut necessary screen only; modular ranks are lower bounds, so possible means not excluded by these cuts"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"family_counts": dict(count_by_size),
                      "rank_profiles": result["modular_rank_profile_counts"],
                      "possible": len(possible)}, sort_keys=True))


if __name__ == "__main__":
    main()
