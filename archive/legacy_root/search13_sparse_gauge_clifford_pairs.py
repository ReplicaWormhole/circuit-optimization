"""One/two sparse complex-Hadamard output rotations over Q(zeta_96).

For any same-eigenvalue row pair (a,b), use the unitary block
[[1,q],[r,-r*q]]/sqrt(2), q,r in {1,-1,i,-i}. Enumerate all one-pair blocks
and all disjoint two-pair blocks, composed on the left of the exact rank-eight
output gauge. Modular ranks at mixed masks 83 and 163 are rigorous lower
bounds; exact characteristic-zero ranks of minimizers are checked separately.
"""

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_fulltail_invariant_exact_ansatz import LABELS, tail_matrix
from search13_fulltail_invariant_exact_witness import gauge_matrix
from search13_rank8_neighborhood_mixedcut import field_mod, rank_mod, realignment
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_sparse_gauge_clifford_pairs_result.json"
MASKS = (83, 163)
PRIME = 97


def pair_groups():
    pairs = [(a, b) for a, b in combinations(range(16), 2)
             if LABELS[a] == LABELS[b]]
    doubles = [group for group in combinations(pairs, 2)
               if len({row for pair in group for row in pair}) == 4]
    assert len(pairs) == 27 and len(doubles) == 273
    return [(pair,) for pair in pairs] + doubles


def target_mod(prime):
    exact = CyclotomicPair()
    g, _ = gauge_matrix(exact)
    target = exact.matmul(g, tail_matrix(exact))
    root = pow(primitive_root(prime), (prime - 1) // 96, prime)
    assert pow(root, 96, prime) == 1 and pow(root, 48, prime) != 1
    matrix = np.asarray([[field_mod(v, root, prime) for v in row]
                         for row in target], dtype=np.int64)
    h = (exact.z**12 + exact.z**-12) * exact.half
    return matrix, field_mod(h, root, prime), pow(root, 24, prime)


def rotate(base, group, phase_indices, h, imag, prime):
    out = base.copy()
    values = (1, prime - 1, imag, (-imag) % prime)
    for (a, b), (q_index, r_index) in zip(group, phase_indices):
        q, r = values[q_index], values[r_index]
        old_a, old_b = out[a].copy(), out[b].copy()
        out[a] = (h * (old_a + q * old_b)) % prime
        out[b] = (h * r * (old_a - q * old_b)) % prime
    return out


def main():
    assert not OUT.exists()
    base, h, imag = target_mod(PRIME)
    assert (h*h*2) % PRIME == 1 and (imag*imag) % PRIME == PRIME - 1
    profiles = Counter()
    minima = {mask: (99, None) for mask in MASKS}
    sum_minimum = (999, None)
    potentially_relevant = []
    potentially_relevant_count = 0
    counts = Counter()
    for group in pair_groups():
        for phase_choice in product(product(range(4), repeat=2), repeat=len(group)):
            transformed = rotate(base, group, phase_choice, h, imag, PRIME)
            ranks = tuple(rank_mod(realignment(transformed, mask), PRIME)
                          for mask in MASKS)
            counts[len(group)] += 1
            profiles[ranks] += 1
            item = {"row_pairs": group, "phase_indices": phase_choice,
                    "mod97_ranks": ranks}
            for mask, rank in zip(MASKS, ranks):
                if rank < minima[mask][0]:
                    minima[mask] = (rank, item)
            if sum(ranks) < sum_minimum[0]:
                sum_minimum = (sum(ranks), item)
            if min(ranks) <= 8:
                potentially_relevant_count += 1
                if len(potentially_relevant) < 100:
                    potentially_relevant.append(item)
    assert counts == {1: 432, 2: 69888}
    result = {"base_gauge": "search13_fulltail_invariant_exact_witness.py",
              "family": "one or two disjoint same-eigenvalue 2x2 complex-Hadamard row blocks [[1,q],[r,-r*q]]/sqrt(2), q,r fourth roots of unity",
              "phase_index_to_value": ["1", "-1", "i", "-i"],
              "field": "Q(zeta_96), rank reduced modulo 97",
              "mask_order": MASKS, "counts": dict(counts),
              "mod97_rank_profiles": {f"{a},{b}": n for (a, b), n in sorted(profiles.items())},
              "minimum_mod97_ranks_by_mask": {str(mask): {"rank": rank, "example": item}
                                               for mask, (rank, item) in minima.items()},
              "minimum_sum_mod97_ranks": {"sum": sum_minimum[0],
                                          "example": sum_minimum[1]},
              "potentially_relevant_rank_at_most_8_count": potentially_relevant_count,
              "up_to_100_potentially_relevant_examples": potentially_relevant,
              "scope": "exact finite-field lower-bound screen for one/two pair rotations at two masks; exact characteristic-zero rank requires separate verification"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"counts": dict(counts),
                      "minimum_mod97_ranks": {str(k): v[0] for k, v in minima.items()},
                      "potentially_relevant": potentially_relevant_count}))


if __name__ == "__main__":
    main()
