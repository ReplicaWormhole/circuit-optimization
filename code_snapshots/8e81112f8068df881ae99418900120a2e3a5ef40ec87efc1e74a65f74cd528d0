"""Exact spectral-path Gram for selected Pauli image on physical chain 0-2-1."""

import json
from itertools import combinations_with_replacement, permutations

import numpy as np

from global_lower_bound6_chain_gram import kron_label, rank_mod_prime


def main():
    labels = [a + c + b + 'I' for a in 'XYZ'
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
    assert int(np.max(np.abs(gram))) < 2 ** 27
    result = {
        'chain': [0, 2, 1], 'gauge': 'final wire 1 axis Z by collective rotation',
        'basis_labels': labels, 'monomials': [list(pair) for pair in monomials],
        'gram': gram.tolist(), 'rank_mod_prime': rank_mod_prime(gram, 1000003),
        'prime': 1000003, 'max_absolute_gram_entry': int(np.max(np.abs(gram))),
    }
    with open('global_lower_bound6_opposite_chain_gram_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('gram', 'monomials', 'basis_labels')}))


if __name__ == '__main__':
    main()
