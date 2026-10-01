"""Exact rank-one obstruction and count for selected physical chain 0-1-3."""

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
from global_lower_bound6_chain_certificate import (
    chain_match, forced_zero_coordinates,
)


def adjacent_v_match(schedule, deg):
    matches = []
    for j in range(4):
        if deg[j] != 2 or j not in schedule[3] or j in schedule[4]:
            continue
        k = next(wire for wire in schedule[3] if wire != j)
        if k not in schedule[4]:
            continue
        ell = next(wire for wire in schedule[4] if wire != k)
        if tuple(sorted((j, k))) in ADJACENT and tuple(sorted((j, ell))) in ADJACENT:
            matches.append((j, k, ell))
    return matches


def main():
    with open('global_lower_bound6_adjacent_v_gram_result.json', encoding='utf-8') as source:
        data = json.load(source)
    labels = data['basis_labels']
    monomials = [tuple(pair) for pair in data['monomials']]
    gram = sp.Matrix(data['gram'])
    assert labels[0] == 'XIII' and labels[-1] == 'ZZIZ'
    assert gram.shape == (231, 231) and gram == gram.T
    assert len(gram.nullspace()) == 21

    forced_zero, stages = forced_zero_coordinates(gram, monomials, len(labels))
    assert [stage['new_forced_zero_axes'] for stage in stages] == [
        [0, 6, 7, 13, 14, 18, 19, 20], [],
    ]
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
    assert len(basis.polys) == 1 and basis.polys[0].as_expr() == 1

    prior = Counter()
    excluded = Counter()
    remaining = Counter()
    examples = []
    for schedule in product(EDGES, repeat=5):
        deg = degrees(schedule)
        if min(deg) < 2 or not full_forward_support(schedule):
            continue
        if any(terminal_obstructions(schedule, deg)):
            continue
        if len(terminal_opposite_degree_two_wires(schedule, deg)) >= 2:
            continue
        if chain_match(schedule, deg):
            continue
        kind = graph_type(schedule, deg)
        prior[kind] += 1
        matches = adjacent_v_match(schedule, deg)
        if matches:
            excluded[kind] += 1
            if len(examples) < 4:
                examples.append({'edges': [list(edge) for edge in schedule],
                                 'selected_chain': list(matches[0])})
        else:
            remaining[kind] += 1
    assert sum(prior.values()) == 208
    assert sum(excluded.values()) == 120
    assert sum(remaining.values()) == 88
    result = {
        'gram_rank_exact': 231 - len(gram.nullspace()),
        'forced_zero_stages': stages,
        'remaining_coefficient_count_before_groebner': len(keep),
        'groebner_basis_with_unit_sphere': ['1'],
        'prior_survivors': dict(sorted(prior.items())),
        'newly_excluded': dict(sorted(excluded.items())),
        'remaining': dict(sorted(remaining.items())),
        'totals': {'prior_survivors': sum(prior.values()),
                   'newly_excluded': sum(excluded.values()),
                   'remaining': sum(remaining.values())},
        'excluded_examples': examples,
    }
    with open('global_lower_bound6_adjacent_v_certificate_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps(result['totals'], sort_keys=True))


if __name__ == '__main__':
    main()
