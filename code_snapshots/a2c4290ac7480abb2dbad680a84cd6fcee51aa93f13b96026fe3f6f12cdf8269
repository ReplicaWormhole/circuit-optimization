"""Exact modular spectral-path screen for pair-supported selected axes.

For the fixed-prefix residual shift A=P V4 P^-1, every tail diagonalizer
image of an input Pauli axis obeys M_a B M_b B M_c=0 for all distinct
spectral labels. The quadratic image map is tested on all 15 traceless
Pauli strings of each physical pair. Full modular column rank certifies
there is no nonzero pair-supported image over characteristic zero.
"""

from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_fulltail_invariant_commutant import (
    pauli_action, prefix_action, residual_shift_action,
)
from global_lower_bound6_chain_gram import rank_mod_prime

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_fulltail_rank8_spectral_pair_result.json"
PRIMES = (97, 193)
PAIRS = tuple(combinations(range(4), 2))


def setup(prime):
    root = pow(primitive_root(prime), (prime - 1) // 16, prime)
    imag = pow(root, 4, prime)
    perm, exponents = residual_shift_action(*prefix_action())
    a = np.zeros((16, 16), dtype=np.int64)
    for x, y in enumerate(perm):
        a[y, x] = pow(root, exponents[x], prime)
    powers = [np.eye(16, dtype=np.int64)]
    for _ in range(3):
        powers.append((a @ powers[-1]) % prime)
    ms = [sum(pow(imag, -k*r, prime) * powers[r]
              for r in range(4)) % prime for k in range(4)]
    return root, imag, ms


def pair_labels(pair):
    labels = ["".join(axes) for axes in product("IXYZ", repeat=4)
              if any(a != "I" for a in axes)
              and all(a == "I" for j, a in enumerate(axes) if j not in pair)]
    assert len(labels) == 15
    return labels


def image_matrix(pair, prime, imaginary, ms):
    labels = pair_labels(pair)
    paulis = []
    for label in labels:
        flip, coeff = pauli_action(label, imaginary, prime)
        matrix = np.zeros((16, 16), dtype=np.int64)
        for x in range(16):
            matrix[x ^ flip, x] = coeff[x]
        paulis.append(matrix)
    monomials = list(combinations_with_replacement(range(15), 2))
    blocks = []
    for a, b, c in permutations(range(4), 3):
        chunks = []
        for i, j in monomials:
            block = ms[a] @ paulis[i] @ ms[b] @ paulis[j] @ ms[c]
            if i != j:
                block += ms[a] @ paulis[j] @ ms[b] @ paulis[i] @ ms[c]
            chunks.append((block % prime).reshape(-1))
        blocks.append(np.stack(chunks))
    return np.concatenate(blocks, axis=1)


def main():
    rows = []
    for prime in PRIMES:
        root, imaginary, ms = setup(prime)
        for pair in PAIRS:
            image = image_matrix(pair, prime, imaginary, ms)
            # A full-rank ordinary bilinear Gram is a sufficient exact
            # witness of full quadratic-image rank. Failure is inconclusive.
            gram = (image @ image.T) % prime
            rank = rank_mod_prime(gram, prime)
            rows.append({"prime": prime, "pair": list(pair),
                         "quadratic_image_count": 120,
                         "ordinary_gram_rank_mod_prime": rank,
                         "full_rank_certificate": rank == 120})
            print(json.dumps(rows[-1]), flush=True)
    result = {"residual_shift": "A=P V4 P^-1 after exact14 six-CNOT prefix",
              "field_primes": PRIMES, "rows": rows,
              "scope": "full modular rank certifies pair no-go; deficient Gram alone does not classify rank-one solutions"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
