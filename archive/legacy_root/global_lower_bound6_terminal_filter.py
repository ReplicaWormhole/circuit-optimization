"""Finite necessary-condition filter for all five-CNOT V4 edge schedules.

This scans only unoriented wire pairs: arbitrary CNOT orientations and local
gates are left free. It proves exclusions only from the previously established
degree, support, and terminal spectral-path lemmas.
"""

from collections import Counter, defaultdict
from itertools import product
import json


EDGES = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
ADJACENT = {(0, 1), (1, 2), (2, 3), (0, 3)}


def degrees(schedule):
    counts = [0] * 4
    for a, b in schedule:
        counts[a] += 1
        counts[b] += 1
    return tuple(counts)


def full_forward_support(schedule):
    """Necessary because U Z_j U† commutes with V4 and has full support."""
    for j in range(4):
        support = {j}
        for a, b in schedule:
            if a in support or b in support:
                support.update((a, b))
        if len(support) != 4:
            return False
    return True


def terminal_obstructions(schedule, deg):
    """Spectral one-CNOT image obstructions for terminal degree-two wires."""
    adjacent_bad = []
    opposite_double_bad = []
    for t, (a, b) in enumerate(schedule):
        if any({a, b} & set(later) for later in schedule[t + 1:]):
            continue
        if deg[a] != 2 and deg[b] != 2:
            continue
        if (a, b) in ADJACENT:
            adjacent_bad.append(t)
        elif deg[a] == deg[b] == 2:
            opposite_double_bad.append(t)
    return adjacent_bad, opposite_double_bad


def main():
    counts = defaultdict(Counter)
    survivors = defaultdict(list)
    for schedule in product(EDGES, repeat=5):
        deg = degrees(schedule)
        if min(deg) < 2:
            continue
        pattern = tuple(sorted(deg, reverse=True))
        assert pattern in ((4, 2, 2, 2), (3, 3, 2, 2))
        label = "".join(map(str, pattern))
        counts[label]["degree_feasible"] += 1
        if not full_forward_support(schedule):
            continue
        counts[label]["full_forward_support"] += 1
        adjacent_bad, opposite_double_bad = terminal_obstructions(schedule, deg)
        if adjacent_bad:
            counts[label]["terminal_adjacent_excluded"] += 1
        if opposite_double_bad:
            counts[label]["terminal_opposite_double_excluded"] += 1
        if adjacent_bad or opposite_double_bad:
            counts[label]["terminal_excluded_union"] += 1
            continue
        counts[label]["surviving_schedule_count"] += 1
        if len(survivors[label]) < 8:
            survivors[label].append({
                "edges": [list(edge) for edge in schedule],
                "degrees": list(deg),
            })
    output = {
        "total_ordered_unoriented_pair_schedules": len(EDGES) ** 5,
        "rules": [
            "Each wire has at least two CNOT incidences.",
            "Forward support of every input Z_j must cover all four physical wires.",
            "Terminal adjacent edge with a degree-two endpoint is impossible.",
            "Terminal opposite edge with two degree-two endpoints is impossible.",
        ],
        "by_degree_pattern": {k: dict(v) for k, v in sorted(counts.items())},
        "survivor_examples": dict(survivors),
    }
    with open("global_lower_bound6_terminal_filter_result.json", "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({k: v for k, v in output.items() if k != "survivor_examples"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
