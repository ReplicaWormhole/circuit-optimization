"""Exact necessary half-shift symmetry for the selected 0-2-1 chain.

Unlike the adjacent-first chain results, the rank-one locus here is
nonempty. The exact ideal nevertheless forces every real spectral-path
solution in the 21D span to be supported on the opposite pair 0,2 and
invariant under V4^2.
"""

from collections import Counter
from itertools import product
import json

import sympy as sp

from global_lower_bound6_terminal_filter import (
    ADJACENT, EDGES, degrees, full_forward_support, terminal_obstructions,
)
from global_lower_bound6_two_terminal_axes import (
    graph_type, terminal_opposite_degree_two_wires,
)
from global_lower_bound6_chain_certificate import chain_match, forced_zero_coordinates
from global_lower_bound6_adjacent_v_certificate import adjacent_v_match


def opposite_chain_match(schedule, deg):
    matches = []
    for j in range(4):
        if deg[j] != 2 or j not in schedule[3] or j in schedule[4]:
            continue
        k = next(wire for wire in schedule[3] if wire != j)
        if k not in schedule[4]:
            continue
        ell = next(wire for wire in schedule[4] if wire != k)
        if tuple(sorted((j, k))) not in ADJACENT:
            matches.append((j, k, ell))
    return matches


def main():
    with open('global_lower_bound6_opposite_chain_gram_result.json', encoding='utf-8') as source:
        data = json.load(source)
    labels = data['basis_labels']
    monomials = [tuple(pair) for pair in data['monomials']]
    gram = sp.Matrix(data['gram'])
    assert labels[0] == 'XIII' and labels[-1] == 'ZZZI'
    assert gram.shape == (231, 231) and gram == gram.T
    assert len(gram.nullspace()) == 42

    forced_zero, stages = forced_zero_coordinates(gram, monomials, len(labels))
    assert sorted(forced_zero) == [0, 6, 7, 13, 14, 18, 19, 20]
    keep = [i for i in range(len(labels)) if i not in forced_zero]
    variables = sp.symbols(f'x0:{len(keep)}')
    variable_at = {original: variables[position] for position, original in enumerate(keep)}
    monomial_vector = sp.Matrix([
        variable_at[i] * variable_at[j] if i in variable_at and j in variable_at else 0
        for i, j in monomials
    ])
    equations = set(gram * monomial_vector)
    equations.discard(sp.Integer(0))
    equations.add(sum(x * x for x in variables) - 1)
    basis = sp.groebner(list(equations), *variables, order='grevlex', domain=sp.QQ)
    three_body = [4, 5, 11, 12]
    antisymmetric_pairs = [(2, 8), (3, 15), (10, 16)]
    for i in three_body:
        assert basis.reduce(variable_at[i] ** 2)[1] == 0
    for i, j in antisymmetric_pairs:
        assert basis.reduce((variable_at[i] - variable_at[j]) ** 2)[1] == 0
    # An exact normalized solution remains: X_0 X_2.
    witness = sp.zeros(len(labels), 1)
    witness[labels.index('XIXI')] = 1
    witness_monomials = sp.Matrix([witness[i] * witness[j] for i, j in monomials])
    assert gram * witness_monomials == sp.zeros(231, 1)

    prior = Counter()
    classified = Counter()
    examples = []
    terminal_axis_counts = Counter()
    for schedule in product(EDGES, repeat=5):
        deg = degrees(schedule)
        if min(deg) < 2 or not full_forward_support(schedule):
            continue
        if any(terminal_obstructions(schedule, deg)):
            continue
        selected = terminal_opposite_degree_two_wires(schedule, deg)
        if len(selected) >= 2 or chain_match(schedule, deg) or adjacent_v_match(schedule, deg):
            continue
        kind = graph_type(schedule, deg)
        prior[kind] += 1
        matches = opposite_chain_match(schedule, deg)
        if matches:
            classified[kind] += 1
            terminal_axis_counts[len(selected)] += 1
            if len(examples) < 4:
                examples.append({'edges': [list(edge) for edge in schedule],
                                 'selected_chain': list(matches[0])})
    assert sum(prior.values()) == 88
    assert sum(classified.values()) == 48
    assert terminal_axis_counts == Counter({0: 48})
    result = {
        'gram_rank_exact': 231 - len(gram.nullspace()),
        'forced_zero_stages': stages,
        'three_body_coefficients_forced_zero': [labels[i] for i in three_body],
        'opposite_pair_coefficient_equalities': [
            [labels[i], labels[j]] for i, j in antisymmetric_pairs],
        'exact_surviving_image': 'X0 X2',
        'prior_survivors': dict(sorted(prior.items())),
        'classified_schedules': dict(sorted(classified.items())),
        'classified_terminal_opposite_axis_counts': dict(sorted(terminal_axis_counts.items())),
        'classified_examples': examples,
    }
    with open('global_lower_bound6_opposite_chain_certificate_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps({'classified': sum(classified.values()),
                      'prior_survivors': sum(prior.values()),
                      'new_exclusions': 0}, sort_keys=True))


if __name__ == '__main__':
    main()
