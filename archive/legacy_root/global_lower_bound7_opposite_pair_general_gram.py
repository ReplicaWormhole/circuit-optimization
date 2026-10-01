"""Exact integer spectral-path Gram on every traceless opposite-pair Pauli."""

from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path

import numpy as np

from global_lower_bound6_chain_gram import kron_label, rank_mod_prime

ROOT = Path(__file__).resolve().parent


def main():
    labels = [a + 'I' + b + 'I' for a in 'IXYZ' for b in 'IXYZ'
              if a + b != 'II']
    assert len(labels) == 15
    axes = [kron_label(label) for label in labels]
    shift = np.zeros((16, 16), dtype=np.complex128)
    for x in range(16):
        shift[((x & 1) << 3) | (x >> 1), x] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in range(4)]
    ms = [sum(((-1j) ** (k*r) * powers[r] for r in range(4)),
              np.zeros((16,16), dtype=np.complex128)) for k in range(4)]
    monomials = list(combinations_with_replacement(range(15), 2))
    blocks = []
    for a,b,c in permutations(range(4),3):
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
    assert int(np.max(np.abs(gram))) < 2**27
    rank = rank_mod_prime(gram, 1000003)
    result = {'pair':'02', 'labels':labels,
              'monomials':[list(x) for x in monomials], 'gram':gram.tolist(),
              'rank_mod_prime':rank, 'prime':1000003,
              'max_absolute_gram_entry':int(np.max(np.abs(gram)))}
    (ROOT / 'global_lower_bound7_opposite_pair_general_gram_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
                      ('labels','monomials','gram')}, sort_keys=True))


if __name__ == '__main__':
    main()
