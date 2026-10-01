"""Depth-seven search after cancelling the first Fourier CZ parity mask.

Adding pi/2 to all four cyclic adjacent-pair Z masks is shift invariant.
It cancels the CZ mask Z1Z2 and creates three other pair masks. The triple
gauge is varied over the four choices that cancel one of its four masks.
This script asks only whether a closed seven-CNOT parity network visits the
six remaining masks. No result here constrains non-diagonalizing bases.
"""

import json
from collections import deque


TRIPLES = (11, 13, 14, 7)
PAIRS = (12, 3, 9)  # 01, 23, 30; 12 cancels


def search_depth_seven(required):
    start = ((8, 4, 2, 1), 0)
    goal = ((8, 4, 2, 1), (1 << len(required)) - 1)
    flags = {mask: 1 << i for i, mask in enumerate(required)}
    queue = deque([(start, 0)])
    previous = {start: None}
    levels = [1] + [0] * 7
    while queue:
        (rows, seen), depth = queue.popleft()
        if (rows, seen) == goal:
            edges = []
            state = goal
            while state != start:
                parent, edge = previous[state]
                edges.append(edge)
                state = parent
            return edges[::-1], levels
        if depth == 7:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                new_rows = list(rows)
                new_rows[target] ^= rows[control]
                new_rows = tuple(new_rows)
                new_seen = seen | flags.get(new_rows[target], 0)
                state = (new_rows, new_seen)
                if state not in previous:
                    previous[state] = ((rows, seen), (control, target))
                    levels[depth + 1] += 1
                    queue.append((state, depth + 1))
    return None, levels


def main():
    results = []
    for cancelled in TRIPLES:
        required = tuple(mask for mask in TRIPLES if mask != cancelled) + PAIRS
        edges, levels = search_depth_seven(required)
        results.append({"cancelled_triple": cancelled,
                        "required_masks": required,
                        "depth_levels": levels,
                        "found_seven_or_less": edges is not None,
                        "edges": edges})
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
