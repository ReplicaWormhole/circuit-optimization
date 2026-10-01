"""Extend run 190's exact five-CNOT schedule filter by a two-axis lemma.

An opposite terminal second incidence of a degree-two wire gives a selected
input Pauli axis whose image commutes with V4^2. Two distinct such wires are
impossible by the (10, 6) multiplicities of V4^2. Only unoriented pairs are
enumerated; all local gates and CNOT directions remain free.
"""

from collections import Counter, defaultdict
from itertools import product
import json

from global_lower_bound6_terminal_filter import (
    ADJACENT, EDGES, degrees, full_forward_support, terminal_obstructions,
)


def terminal_opposite_degree_two_wires(schedule, deg):
    selected = []
    for wire in range(4):
        if deg[wire] != 2:
            continue
        second = [t for t, edge in enumerate(schedule) if wire in edge][1]
        edge = schedule[second]
        if edge in ADJACENT:
            continue
        if any(set(edge) & set(later) for later in schedule[second + 1:]):
            continue
        selected.append(wire)
    return tuple(selected)


def graph_type(schedule, deg):
    pattern = ''.join(map(str, sorted(deg, reverse=True)))
    multiplicities = ''.join(map(str, sorted(Counter(schedule).values())))
    return f'{pattern}:{multiplicities}'


def main():
    counts = Counter()
    eliminated = Counter()
    survivors = Counter()
    examples = defaultdict(list)
    for schedule in product(EDGES, repeat=5):
        deg = degrees(schedule)
        if min(deg) < 2 or not full_forward_support(schedule):
            continue
        if any(terminal_obstructions(schedule, deg)):
            continue
        kind = graph_type(schedule, deg)
        counts[kind] += 1
        selected = terminal_opposite_degree_two_wires(schedule, deg)
        if len(selected) >= 2:
            eliminated[kind] += 1
            if len(examples[kind]) < 2:
                examples[kind].append({
                    'edges': [list(e) for e in schedule],
                    'selected_input_wires': list(selected),
                })
        else:
            survivors[kind] += 1
    assert sum(counts.values()) == 344
    assert sum(eliminated.values()) == 88
    assert sum(survivors.values()) == 256
    result = {
        'prior_survivors': dict(sorted(counts.items())),
        'newly_excluded': dict(sorted(eliminated.items())),
        'remaining': dict(sorted(survivors.items())),
        'totals': {
            'prior_survivors': sum(counts.values()),
            'newly_excluded': sum(eliminated.values()),
            'remaining': sum(survivors.values()),
        },
        'excluded_examples': dict(examples),
    }
    with open('global_lower_bound6_two_terminal_axes_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps(result['totals'], sort_keys=True))


if __name__ == '__main__':
    main()
