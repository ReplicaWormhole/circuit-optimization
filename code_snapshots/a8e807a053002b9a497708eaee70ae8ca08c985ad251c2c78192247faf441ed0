"""Exact necessary-condition screen for all six-CNOT pair schedules.

This uses only obstructions proved in global_lower_bound6_complete_filter.md
and its referenced certificates.  A survivor is not a valid diagonalizer.
"""

from collections import Counter
from itertools import product
import json
from pathlib import Path

EDGES = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
OPPOSITE = {(0, 2), (1, 3)}
ROOT = Path(__file__).resolve().parent


def degrees(schedule):
    result = [0] * 4
    for a, b in schedule:
        result[a] += 1
        result[b] += 1
    return tuple(result)


def all_forward_cones_full(schedule):
    for j in range(4):
        cone = {j}
        for edge in schedule:
            if cone.intersection(edge):
                cone.update(edge)
        if len(cone) != 4:
            return False
    return True


def selected_tail(schedule, j):
    """Potential support after the second incidence of degree-two wire j.

    The selected input axis commutes through its first incidence.  CNOTs
    disjoint from its potential support then act trivially on it.  This
    deliberately overapproximates support, independently of orientations.
    """
    second = [t for t, edge in enumerate(schedule) if j in edge][1]
    first_edge = schedule[second]
    cone = set(first_edge)
    active = []
    for t in range(second + 1, len(schedule)):
        edge = schedule[t]
        if cone.intersection(edge):
            active.append((t, edge))
            cone.update(edge)
    return second, first_edge, active


def tail_rejection(schedule, deg):
    """Return one proven exclusion reason, or None.

    A zero-active-gate tail is the terminal-pair theorem.  With exactly one
    active gate, the selected image belongs to the prior 21-dimensional
    two-CNOT-chain overapproximation.  Any disjoint intervening gates act
    as the identity on that selected image.
    """
    for j in range(4):
        if deg[j] != 2:
            continue
        _, edge, active = selected_tail(schedule, j)
        if not active:
            return ('terminal_opposite_trace' if edge in OPPOSITE
                    else 'terminal_adjacent_gram'), j
        if len(active) != 1:
            continue
        _, later = active[0]
        k = next(w for w in edge if w != j)
        # j has degree two, so the only possible active endpoint is k.
        assert k in later and j not in later
        ell = next(w for w in later if w != k)
        if edge in OPPOSITE:
            return 'one_active_opposite_first_trace', j
        if tuple(sorted((j, ell))) in OPPOSITE:
            return 'one_active_adjacent_chain', j
        return 'one_active_adjacent_v', j
    return None


def classify(length):
    counts = Counter()
    survivor_degrees = Counter()
    survivor_tail_profiles = Counter()
    examples = {}
    for schedule in product(EDGES, repeat=length):
        deg = degrees(schedule)
        if min(deg) < 2:
            reason = 'degree_rejected'
        elif not all_forward_cones_full(schedule):
            reason = 'forward_support_rejected'
        else:
            rejection = tail_rejection(schedule, deg)
            reason = rejection[0] if rejection else 'survives_known_tests'
        counts[reason] += 1
        if reason == 'survives_known_tests':
            survivor_degrees[''.join(map(str, sorted(deg, reverse=True)))] += 1
            profile = tuple(sorted(len(selected_tail(schedule, j)[2])
                                   for j in range(4) if deg[j] == 2))
            survivor_tail_profiles[','.join(map(str, profile)) or 'none'] += 1
        examples.setdefault(reason, [list(edge) for edge in schedule])
    assert sum(counts.values()) == 6 ** length
    return {
        'cnot_count': length,
        'ordered_unoriented_pair_schedules': 6 ** length,
        'exclusive_partition': dict(sorted(counts.items())),
        'survivor_degree_sequences': dict(sorted(survivor_degrees.items())),
        'survivor_active_tail_profiles': dict(sorted(survivor_tail_profiles.items())),
        'examples': examples,
    }


def main():
    five = classify(5)
    assert five['exclusive_partition'].get('survives_known_tests', 0) == 0
    six = classify(6)
    result = {'five_cnot_regression': five, 'six_cnot_screen': six}
    path = ROOT / 'global_lower_bound7_schedule_filter_result.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'five': five['exclusive_partition'],
                      'six': six['exclusive_partition']}, sort_keys=True))


if __name__ == '__main__':
    main()
