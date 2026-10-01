"""Exact spectral-path Gram for a selected axis on physical chain 0-1-2.

The axis commutes through its first incident CNOT, meets CNOT(0,1) at its
second incidence, and then meets CNOT(1,2), with no subsequent CNOT. A
collective rotation puts the final factor on wire 2 along Z. Every resulting
image lies in the real span of the 21 Pauli strings listed below.

Entries are exact Gaussian integers despite complex128 arithmetic: each
projector numerator has Gaussian-integer entries; each quadratic spectral
image entry has magnitude <=128, and each Gram entry has magnitude <2^27.
"""

import json
from itertools import combinations_with_replacement, permutations

import numpy as np


PAULI = {
    'I': np.array([[1, 0], [0, 1]], dtype=np.complex128),
    'X': np.array([[0, 1], [1, 0]], dtype=np.complex128),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    'Z': np.array([[1, 0], [0, -1]], dtype=np.complex128),
}


def kron_label(label):
    out = np.array([[1]], dtype=np.complex128)
    for letter in label:
        out = np.kron(out, PAULI[letter])
    return out


def rank_mod_prime(matrix, prime):
    a = np.mod(matrix, prime).astype(np.int64)
    n = len(a)
    rank = 0
    for col in range(n):
        pivot = next((row for row in range(rank, n) if a[row, col]), None)
        if pivot is None:
            continue
        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, prime) % prime
        for row in range(rank + 1, n):
            if a[row, col]:
                a[row] = (a[row] - a[row, col] * a[rank]) % prime
        rank += 1
        if rank == n:
            break
    return rank


def main():
    labels = [a + b + c + 'I' for a in 'XYZ'
              for b, c in [('I', 'I'), *[(b, 'I') for b in 'XYZ'],
                           *[(b, 'Z') for b in 'XYZ']]]
    assert len(labels) == 21 and len(set(labels)) == 21
    axes = [kron_label(label) for label in labels]
    shift = np.zeros((16, 16), dtype=np.complex128)
    for x in range(16):
        shift[((x & 1) << 3) | (x >> 1), x] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in range(4)]
    projectors4 = [sum(((-1j) ** (k * r) * powers[r] for r in range(4)),
                       np.zeros((16, 16), dtype=np.complex128)) for k in range(4)]
    monomials = list(combinations_with_replacement(range(len(axes)), 2))
    rows = []
    for a, b, c in permutations(range(4), 3):
        for i, j in monomials:
            block = projectors4[a] @ axes[i] @ projectors4[b] @ axes[j] @ projectors4[c]
            if i != j:
                block += projectors4[a] @ axes[j] @ projectors4[b] @ axes[i] @ projectors4[c]
            rows.append(block.reshape(-1))
    image = np.array(rows).reshape(24, len(monomials), 256).transpose(1, 0, 2).reshape(len(monomials), -1)
    gram_complex = image.conj() @ image.T
    assert np.max(np.abs(gram_complex.imag)) == 0
    assert np.max(np.abs(gram_complex.real - np.rint(gram_complex.real))) == 0
    gram = np.rint(gram_complex.real).astype(np.int64)
    assert np.array_equal(gram, gram.T)
    max_entry = int(np.max(np.abs(gram)))
    assert max_entry < 2 ** 27
    prime = 1000003
    rank = rank_mod_prime(gram, prime)
    result = {
        'chain': [0, 1, 2], 'gauge': 'final wire 2 axis Z by collective rotation',
        'basis_labels': labels, 'monomials': [list(pair) for pair in monomials],
        'gram': gram.tolist(), 'rank_mod_prime': rank, 'prime': prime,
        'max_absolute_gram_entry': max_entry,
    }
    with open('global_lower_bound6_chain_gram_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('gram', 'monomials', 'basis_labels')}))


if __name__ == '__main__':
    main()
