"""Exclude schedules with a first-commuting selected axis confined to an adjacent pair.

Uses the exact arbitrary-adjacent-pair spectral-path certificate from run 240.
Potential support is an orientation-independent overapproximation.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

import global_lower_bound7_schedule_filter as prior

ROOT = Path(__file__).resolve().parent


def first_commuting_cone(schedule, j):
    first = next(t for t, edge in enumerate(schedule) if j in edge)
    cone = {j}
    active = []
    for t in range(first + 1, len(schedule)):
        edge = schedule[t]
        if cone.intersection(edge):
            cone.update(edge)
            active.append(t)
    return cone, active


def main():
    counts = Counter()
    excluded_by_degree = Counter()
    remaining_by_degree = Counter()
    excluded_examples = {}
    for schedule in product(prior.EDGES, repeat=6):
        deg = prior.degrees(schedule)
        if min(deg) < 2 or not prior.all_forward_cones_full(schedule):
            continue
        if prior.tail_rejection(schedule, deg) is not None:
            continue
        degree_key = ''.join(map(str, sorted(deg, reverse=True)))
        witnesses = []
        for j in range(4):
            cone, active = first_commuting_cone(schedule, j)
            if len(cone) == 2 and tuple(sorted(cone)) not in prior.OPPOSITE:
                witnesses.append({'wire': j, 'pair': sorted(cone),
                                  'active_gate_indices_zero_based': active})
        if witnesses:
            counts['excluded_adjacent_pair'] += 1
            excluded_by_degree[degree_key] += 1
            excluded_examples.setdefault(degree_key, {
                'schedule': [list(e) for e in schedule],
                'selected_axis_witnesses': witnesses})
        else:
            counts['survives_all_applied_tests'] += 1
            remaining_by_degree[degree_key] += 1
    assert sum(counts.values()) == 3576
    result = {
        'prior_survivors': 3576,
        'exclusive_partition': dict(sorted(counts.items())),
        'excluded_by_degree_sequence': dict(sorted(excluded_by_degree.items())),
        'remaining_by_degree_sequence': dict(sorted(remaining_by_degree.items())),
        'excluded_examples': excluded_examples,
    }
    (ROOT / 'global_lower_bound7_adjacent_pair_support_filter_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:result[k] for k in ('exclusive_partition',
                      'excluded_by_degree_sequence','remaining_by_degree_sequence')},
                     sort_keys=True))


if __name__ == '__main__':
    main()
