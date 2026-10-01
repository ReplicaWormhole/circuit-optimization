"""Exact Gram test of the spectral-path map on adjacent two-wire axes.

By collective SU(2) symmetry, the partner factor of a one-CNOT image can be
put on Z. The six-dimensional real span below then contains every such image
on physical wires 0,1. This computes a rational integer Gram matrix for the
quadratic spectral-path map, with projectors represented as Mk=4Pk.
"""

import json
from itertools import combinations_with_replacement, permutations

import numpy as np
import sympy as sp


def kron_all(parts):
    out = np.array([[1]], dtype=np.complex128)
    for part in parts:
        out = np.kron(out, part)
    return out


def main():
    I = np.array([[1, 0], [0, 1]], dtype=np.complex128)
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    axes = [kron_all((p, q, I, I)) for q in (I, Z) for p in (X, Y, Z)]
    shift = np.zeros((16, 16), dtype=np.complex128)
    for x in range(16):
        shift[((x & 1) << 3) | (x >> 1), x] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in range(4)]
    projector4 = [sum(((-1j) ** (k * r) * powers[r] for r in range(4)),
                      np.zeros((16, 16), dtype=np.complex128)) for k in range(4)]
    triples = list(permutations(range(4), 3))
    monomials = list(combinations_with_replacement(range(6), 2))
    images = []
    for a, b, c in triples:
        per_triple = []
        for i, j in monomials:
            block = projector4[a] @ axes[i] @ projector4[b] @ axes[j] @ projector4[c]
            if i != j:
                block += projector4[a] @ axes[j] @ projector4[b] @ axes[i] @ projector4[c]
            per_triple.append(block.reshape(-1))
        images.append(np.stack(per_triple))
    images = np.concatenate(images, axis=1)
    gram_float = (images.conj() @ images.T).real
    assert np.max(np.abs(gram_float - np.rint(gram_float))) == 0
    gram = sp.Matrix(np.rint(gram_float).astype(np.int64).tolist())
    exact_rank = gram.rank()
    eig_min = float(np.linalg.eigvalsh(gram_float)[0])
    result = {"pair": "01", "basis": ["X0", "Y0", "Z0", "X0Z1", "Y0Z1", "Z0Z1"],
              "quadratic_monomial_count": len(monomials), "gram_rank": exact_rank,
              "minimum_numeric_eigenvalue": eig_min,
              "gram": [list(map(int, gram.row(i))) for i in range(gram.rows)],
              "monomials": [list(m) for m in monomials]}
    with open("spectral_path4_adjacent_gram_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("gram", "monomials")}))


if __name__ == "__main__":
    main()
