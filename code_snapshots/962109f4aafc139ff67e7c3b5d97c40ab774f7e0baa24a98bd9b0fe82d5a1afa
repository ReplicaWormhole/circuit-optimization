"""Exact topology inventory for degree-three selected axes in LB7 survivors.

This script applies only the existing one- and two-CNOT certificates when
their hypotheses literally hold.  Other tail shapes are classified, not
excluded.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

import global_lower_bound7_schedule_filter as prior

ROOT = Path(__file__).resolve().parent


def first_commuting_tail(schedule, j):
    first = next(t for t, edge in enumerate(schedule) if j in edge)
    support = {j}
    active = []
    for t in range(first + 1, len(schedule)):
        edge = schedule[t]
        if support.intersection(edge):
            active.append((t, edge))
            support.update(edge)
    return active, support


def geometry(j, active):
    if len(active) != 2:
        return None
    first, second = (edge for _, edge in active)
    assert j in first and j in second
    k = next(w for w in first if w != j)
    ell = next(w for w in second if w != j)
    if k == ell:
        return ('repeated_opposite_pair' if tuple(sorted((j, k))) in prior.OPPOSITE
                else 'repeated_adjacent_pair')
    if tuple(sorted((k, ell))) in prior.OPPOSITE:
        return 'j_center_adjacent_v'
    return 'j_center_mixed_adjacent_opposite'


def main():
    degree3_schedules = 0
    per_wire = Counter()
    schedule_best = Counter()
    degree3333_best = Counter()
    examples = {}
    for schedule in product(prior.EDGES, repeat=6):
        deg = prior.degrees(schedule)
        if min(deg) < 2 or not prior.all_forward_cones_full(schedule):
            continue
        if prior.tail_rejection(schedule, deg) is not None:
            continue
        if any(d == 3 for d in deg):
            degree3_schedules += 1
        tails = []
        for j in range(4):
            if deg[j] != 3:
                continue
            active, support = first_commuting_tail(schedule, j)
            count = len(active)
            if count <= 2:
                kind = geometry(j, active) if count == 2 else f'{count}_active'
            else:
                kind = f'{count}_active'
            key = f'{count}:{kind}'
            tails.append((count, key))
            per_wire[key] += 1
            examples.setdefault(key, {'schedule': [list(e) for e in schedule],
                                      'selected_wire': j,
                                      'active_edges': [list(e) for _, e in active],
                                      'potential_support': sorted(support)})
        if tails:
            best = min(tails)[1]
            schedule_best[best] += 1
            if tuple(sorted(deg)) == (3, 3, 3, 3):
                degree3333_best[best] += 1
    assert degree3_schedules == 3432
    assert sum(degree3333_best.values()) == 1416
    result = {
        'prior_survivors_with_degree_three_wire': degree3_schedules,
        'degree_three_selected_wires_by_tail': dict(sorted(per_wire.items())),
        'surviving_schedules_by_shortest_degree_three_tail': dict(sorted(schedule_best.items())),
        'degree_3333_schedules_by_shortest_tail': dict(sorted(degree3333_best.items())),
        'examples': examples,
        'interpretation': 'classification only; no new schedules excluded',
    }
    (ROOT / 'global_lower_bound7_degree3_tail_screen_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: result[key] for key in
                      ('degree_three_selected_wires_by_tail',
                       'degree_3333_schedules_by_shortest_tail')}, sort_keys=True))


if __name__ == '__main__':
    main()
