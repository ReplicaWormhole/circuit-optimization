"""Exact spectral-path Grams for degree-three two-active j-centered tails.

After an input axis commutes through its first incident CNOT, its only two
active later CNOTs join (j,k) and (j,l), in that order.  The final image
belongs to the 24-dimensional span with Pauli factors on j, optional
arbitrary Pauli on k, and optional fixed-axis Pauli on l.  A collective
rotation sets the last axis on l to Z.  This is an overapproximation of
all local-gate choices and CNOT directions.
"""

from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path

import numpy as np

from global_lower_bound6_chain_gram import kron_label, rank_mod_prime

ROOT = Path(__file__).resolve().parent
SHAPES = {
    'adjacent_adjacent': (0, 1, 3),
    'adjacent_opposite': (0, 1, 2),
    'opposite_adjacent': (0, 2, 1),
}


def basis_labels(j, k, ell):
    labels = []
    for a in 'XYZ':
        for b in 'IXYZ':
            for c in 'IZ':
                label = ['I'] * 4
                label[j], label[k], label[ell] = a, b, c
                labels.append(''.join(label))
    assert len(labels) == 24 and len(set(labels)) == 24
    return labels


def projector_numerators():
    shift = np.zeros((16, 16), dtype=np.complex128)
    for state in range(16):
        shift[((state & 1) << 3) | (state >> 1), state] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in range(4)]
    return [sum(((-1j) ** (a*r) * powers[r] for r in range(4)),
                np.zeros((16,16), dtype=np.complex128)) for a in range(4)]


def gram_for_shape(shape, triple, ms):
    labels = basis_labels(*triple)
    axes = [kron_label(label) for label in labels]
    monomials = list(combinations_with_replacement(range(24), 2))
    blocks = []
    for a,b,c in permutations(range(4), 3):
        per_triple = []
        for i,j in monomials:
            block = ms[a] @ axes[i] @ ms[b] @ axes[j] @ ms[c]
            if i != j:
                block += ms[a] @ axes[j] @ ms[b] @ axes[i] @ ms[c]
            per_triple.append(block.reshape(-1))
        blocks.append(np.stack(per_triple))
    image = np.concatenate(blocks, axis=1)
    gram_float = image.conj() @ image.T
    assert np.max(np.abs(gram_float.imag)) == 0
    assert np.max(np.abs(gram_float.real - np.rint(gram_float.real))) == 0
    gram = np.rint(gram_float.real).astype(np.int64)
    assert np.array_equal(gram, gram.T)
    max_entry = int(np.max(np.abs(gram)))
    assert max_entry < 2**27
    result = {
        'shape': shape, 'ordered_triple_j_k_l': list(triple),
        'basis_labels': labels, 'monomials': [list(pair) for pair in monomials],
        'gram': gram.tolist(), 'rank_mod_prime_1000003': rank_mod_prime(gram, 1000003),
        'quadratic_monomial_count': len(monomials),
        'max_absolute_gram_entry': max_entry,
    }
    path = ROOT / f'global_lower_bound7_degree3_centered_gram_{shape}_result.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return {k:v for k,v in result.items() if k not in ('basis_labels','monomials','gram')}


def main():
    ms = projector_numerators()
    result = {shape:gram_for_shape(shape,triple,ms)
              for shape,triple in SHAPES.items()}
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
