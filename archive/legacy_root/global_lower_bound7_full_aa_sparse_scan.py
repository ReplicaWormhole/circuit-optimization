"""Bounded exact sparse scan of the AA 24D first-axis overapproximation.

Search one-, two-, and three-term signed equal-weight Pauli involutions.
The result classifies only this finite sparse ansatz, not the full 24D span.
"""

from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from global_lower_bound6_chain_gram import kron_label

ROOT = Path(__file__).resolve().parent


def anticommutes(a, b):
    return sum(x != y for x, y in zip(a, b) if x != 'I' and y != 'I') % 2 == 1


def exact_complex_pair(z):
    assert z.real.is_integer() and z.imag.is_integer()
    return int(z.real), int(z.imag)


def main():
    data = json.loads((ROOT / 'global_lower_bound7_degree3_centered_gram_adjacent_adjacent_result.json').read_text())
    labels = data['basis_labels']
    monomials = [tuple(x) for x in data['monomials']]
    gram = np.array(data['gram'], dtype=np.int64)
    assert len(labels) == 24 and gram.shape == (300, 300)
    monomial_index = {m: i for i, m in enumerate(monomials)}
    shift = np.zeros((16, 16), dtype=np.complex128)
    for state in range(16):
        shift[((state & 1) << 3) | (state >> 1), state] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in (1, 2, 3)]
    paulis = [kron_label(label) for label in labels]
    trace_rows = [[exact_complex_pair(np.trace(power @ pauli)) for pauli in paulis]
                  for power in powers]
    counts = {'1': 0, '2': 0, '3': 0}
    solutions = []
    for size in (1, 2, 3):
        for support in combinations(range(24), size):
            if not all(anticommutes(labels[i], labels[j])
                       for i, j in combinations(support, 2)):
                continue
            # Fix global sign, since p and -p define the same input axis.
            for tail in product((-1, 1), repeat=size - 1):
                signs = (1,) + tail
                active = [(monomial_index[(i, j)], signs[a] * signs[b])
                          for a, i in enumerate(support)
                          for b, j in enumerate(support[a:], start=a)]
                norm = sum(ci * cj * int(gram[i, j])
                           for i, ci in active for j, cj in active)
                if norm != 0:
                    continue
                counts[str(size)] += 1
                moments = [tuple(sum(signs[k] * row[i][part]
                                     for k, i in enumerate(support))
                                 for part in (0, 1)) for row in trace_rows]
                invariants = [sum(z*z for z in pair) for pair in moments]
                solutions.append({'labels': [labels[i] for i in support],
                                  'signs': list(signs),
                                  'unnormalized_moments_real_imag': [list(x) for x in moments],
                                  'normalized_moment_abs_squared': [f'{n}/{size}' for n in invariants]})
    canonical = next(s for s in solutions
                     if s['labels'] == ['XIIZ', 'YXIZ', 'ZXII']
                     and s['signs'] == [1, 1, 1])
    canonical_invariant = canonical['normalized_moment_abs_squared']
    distinct = next((s for s in solutions
                     if s['normalized_moment_abs_squared'] != canonical_invariant), None)
    result = {'ansatz': 'at most three pairwise anticommuting Pauli terms, equal absolute coefficients',
              'counts_by_size': counts, 'total_solutions_mod_global_sign': len(solutions),
              'canonical': canonical, 'different_trace_invariant_example': distinct,
              'distinct_trace_invariants': sorted({tuple(s['normalized_moment_abs_squared'])
                                                   for s in solutions}),
              'solutions': solutions}
    (ROOT / 'global_lower_bound7_full_aa_sparse_scan_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'solutions'}, sort_keys=True))


if __name__ == '__main__':
    main()
