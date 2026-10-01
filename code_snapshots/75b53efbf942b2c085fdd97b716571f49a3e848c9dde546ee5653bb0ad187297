"""Bounded exact four-term AA selected-axis support classification.

The 24D adjacent/adjacent span is a sound overapproximation. Enumerate
all four-string supports, retain only supports whose Pauli cross terms
can cancel in A^2, then solve spectral-path and involution equations
with exact modular rank and rational Groebner calculations. The run
limit is explicit; no exclusion is claimed for unprocessed supports.
"""

from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import numpy as np
import sympy as sp

from global_lower_bound6_chain_gram import kron_label, rank_mod_prime
from global_lower_bound7_joint_aa_centralizer import product as pauli_product

ROOT = Path(__file__).resolve().parent
MAX_SINGULAR_SUPPORTS = 250


def cross_groups(strings):
    groups = defaultdict(list)
    for i, j in combinations(range(4), 2):
        phase, label = pauli_product(strings[i], strings[j])
        if phase in (1, -1):
            assert label != 'IIII'
            groups[label].append((i, j, int(phase)))
    return groups


def main():
    source = ROOT / 'global_lower_bound7_degree3_centered_gram_adjacent_adjacent_result.json'
    data = json.loads(source.read_text())
    labels = data['basis_labels']
    monomials = [tuple(pair) for pair in data['monomials']]
    gram = np.asarray(data['gram'], dtype=np.int64)
    assert len(labels) == 24 and gram.shape == (300, 300)
    index = {pair: i for i, pair in enumerate(monomials)}
    shift = np.zeros((16, 16), dtype=np.complex128)
    for state in range(16):
        shift[((state & 1) << 3) | (state >> 1), state] = 1
    paulis = [kron_label(label) for label in labels]
    trace_rows = []
    for r in (1, 2, 3):
        power = np.linalg.matrix_power(shift, r)
        row = []
        for pauli in paulis:
            value = np.trace(power @ pauli)
            assert value.real.is_integer() and value.imag.is_integer()
            row.append((int(value.real), int(value.imag)))
        trace_rows.append(row)

    variables = sp.symbols('u1 u2 u3')
    x = (sp.Integer(1),) + variables
    local_pairs = list(combinations_with_replacement(range(4), 2))
    q = sp.Matrix([x[i]*x[j] for i, j in local_pairs])
    counts = Counter()
    nonfull = []
    complete = True
    for support in combinations(range(24), 4):
        counts['all_four_string_supports'] += 1
        strings = [labels[i] for i in support]
        groups = cross_groups(strings)
        if any(len(terms) == 1 for terms in groups.values()):
            counts['unique_uncancellable_cross_term'] += 1
            continue
        assert all(len(terms) == 2 for terms in groups.values())
        counts['involution_combinatorially_possible'] += 1
        active = [index[(support[i], support[j])] for i, j in local_pairs]
        sub = gram[np.ix_(active, active)]
        if rank_mod_prime(sub, 1000003) == 10:
            counts['spectral_full_rank_mod_prime'] += 1
            continue
        counts['spectral_modular_singular'] += 1
        if counts['spectral_modular_singular'] > MAX_SINGULAR_SUPPORTS:
            complete = False
            break
        matrix = sp.Matrix(sub.tolist())
        nullity = 10 - matrix.rank()
        counts[f'exact_nullity_{nullity}'] += 1
        assert nullity > 0
        involution = [sum(phase*x[i]*x[j] for i, j, phase in terms)
                      for terms in groups.values()]
        equations = list(matrix*q) + involution
        gb = sp.groebner(equations, *variables, order='lex', domain=sp.QQ)
        if gb.contains(sp.Integer(1)):
            counts['joint_ideal_unit'] += 1
            continue
        witness = sp.symbols('w')
        saturated = sp.groebner(equations +
                                [witness*variables[0]*variables[1]*variables[2]-1],
                                *variables, witness, order='lex', domain=sp.QQ)
        if saturated.contains(sp.Integer(1)):
            counts['only_zero_coefficient_solutions'] += 1
            continue
        counts['complex_feasible_nonzero_supports'] += 1
        nonfull.append({'labels': strings, 'exact_gram_nullity': nullity,
                        'commuting_cross_groups':
                        {name: [list(term) for term in terms]
                         for name, terms in sorted(groups.items())},
                        'projective_groebner_basis':
                        [str(p.as_expr()) for p in gb.polys],
                        'nonzero_saturated_groebner_basis':
                        [str(p.as_expr()) for p in saturated.polys],
                        'trace_rows_real_imag':
                        [[list(row[i]) for i in support] for row in trace_rows]})
    result = {'complete_support_enumeration': complete,
              'max_modular_singular_supports_to_solve': MAX_SINGULAR_SUPPORTS,
              'counts': dict(sorted(counts.items())),
              'complex_feasible_nonzero_supports': nonfull,
              'scope': 'all exactly-four nonzero real-Pauli-coefficient supports in fixed-gauge 24D AA span; complex feasibility is not real feasibility'}
    (ROOT / 'global_lower_bound7_full_aa_four_term_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'complete_support_enumeration': complete,
                      'counts': result['counts'],
                      'complex_feasible_nonzero_support_count': len(nonfull)},
                     sort_keys=True))


if __name__ == '__main__':
    main()
