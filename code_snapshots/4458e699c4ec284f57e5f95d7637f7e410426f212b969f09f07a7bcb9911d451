"""Exact bounded classification of all <=3-term AA first-axis involutions.

For at most three distinct Pauli terms, an involution requires pairwise
anticommutation. The spectral-path Gram then gives quadratic equations
in the real coefficients. Normalize the first nonzero coefficient to 1
and solve each support's projective equations by rational Groebner basis.
"""

from collections import Counter
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import numpy as np
import sympy as sp

from global_lower_bound6_chain_gram import rank_mod_prime

ROOT = Path(__file__).resolve().parent


def anticommutes(a, b):
    return sum(x != y for x, y in zip(a, b) if x != 'I' and y != 'I') % 2 == 1


def main():
    data = json.loads((ROOT / 'global_lower_bound7_degree3_centered_gram_adjacent_adjacent_result.json').read_text())
    labels = data['basis_labels']
    gram = np.asarray(data['gram'], dtype=np.int64)
    monomials = [tuple(pair) for pair in data['monomials']]
    assert labels and len(labels) == 24 and gram.shape == (300, 300)
    index = {pair: i for i, pair in enumerate(monomials)}
    counts = Counter()
    feasible = []
    for size in (1, 2, 3):
        for support in combinations(range(24), size):
            if not all(anticommutes(labels[i], labels[j])
                       for i, j in combinations(support, 2)):
                continue
            counts[f'{size}_anticommuting_supports'] += 1
            local_pairs = list(combinations_with_replacement(range(size), 2))
            active = [index[(support[a], support[b])] for a, b in local_pairs]
            sub = gram[np.ix_(active, active)]
            if rank_mod_prime(sub, 1000003) == len(active):
                counts[f'{size}_full_rank_mod_prime'] += 1
                continue
            matrix = sp.Matrix(sub.tolist())
            nullity = len(active) - matrix.rank()
            counts[f'{size}_singular_exact'] += 1
            counts[f'{size}_nullity_{nullity}'] += 1
            variables = sp.symbols('u:'+str(size - 1)) if size > 1 else ()
            x = (sp.Integer(1),) + variables
            q = sp.Matrix([x[a] * x[b] for a, b in local_pairs])
            equations = [sp.expand(z) for z in matrix*q if z != 0]
            if size == 1:
                assert not equations
                basis = []
                has_complex_solution = True
            else:
                gb = sp.groebner(equations, *variables, order='lex', domain=sp.QQ)
                basis = [str(poly.as_expr()) for poly in gb.polys]
                has_complex_solution = not gb.contains(sp.Integer(1))
            if has_complex_solution:
                counts[f'{size}_complex_feasible_supports'] += 1
                feasible.append({'labels': [labels[i] for i in support],
                                 'nullity': nullity,
                                 'projective_groebner_basis': basis})
            else:
                counts[f'{size}_singular_but_infeasible'] += 1
    result = {'scope': 'all <=3 distinct nonzero Pauli coefficients, pairwise anticommuting, in 24D AA span',
              'normalization': 'first nonzero coefficient fixed to 1; real solutions need later real-root check',
              'counts': dict(sorted(counts.items())),
              'complex_feasible_supports': feasible}
    (ROOT / 'global_lower_bound7_full_aa_three_term_real_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
