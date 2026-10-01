"""Apply the exact i-sector Schmidt crossing bound to six-CNOT schedules.

The cut theorem is derived in global_lower_bound7_cut_schmidt.md.  Counts
are applied after the earlier orientation-independent selected-tail filter.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

import global_lower_bound7_schedule_filter as prior

ROOT = Path(__file__).resolve().parent
ADJACENT_CUTS = ((0, 1), (0, 3))


def crossing_count(schedule, left):
    half = set(left)
    return sum((a in half) != (b in half) for a, b in schedule)


def prior_survivor(schedule):
    deg = prior.degrees(schedule)
    return (min(deg) >= 2 and prior.all_forward_cones_full(schedule)
            and prior.tail_rejection(schedule, deg) is None)


def main():
    partition = Counter()
    by_degrees = Counter()
    examples = {}
    for schedule in product(prior.EDGES, repeat=6):
        if not prior_survivor(schedule):
            continue
        n01, n03 = (crossing_count(schedule, left) for left in ADJACENT_CUTS)
        if n01 < 2 and n03 < 2:
            category = 'both_adjacent_cuts_below_two'
        elif n01 < 2:
            category = 'cut_01_below_two'
        elif n03 < 2:
            category = 'cut_03_below_two'
        else:
            category = 'survives_cut_and_tail_tests'
        partition[category] += 1
        degree_key = ''.join(map(str, sorted(prior.degrees(schedule), reverse=True)))
        by_degrees[(category, degree_key)] += 1
        examples.setdefault(category, [list(edge) for edge in schedule])
    assert sum(partition.values()) == 3576
    result = {
        'prior_six_cnot_survivors': 3576,
        'exclusive_partition': dict(sorted(partition.items())),
        'by_degree_sequence': {
            category: {degree: by_degrees[(category, degree)]
                       for degree in ('3333', '4332', '4422', '5322')
                       if by_degrees[(category, degree)]}
            for category in sorted(partition)
        },
        'examples': examples,
    }
    path = ROOT / 'global_lower_bound7_cut_crossing_filter_result.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result['exclusive_partition'], sort_keys=True))


if __name__ == '__main__':
    main()
