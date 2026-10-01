"""Enumerate six-CNOT parity prefixes with a two-CNOT endpoint map.

This is a structural search within the parity-rotation formulation, not a
lower bound for arbitrary 15-CNOT diagonalizers. It records which of the four
required parity masks have appeared on a target wire after each CNOT.
"""

import json
from collections import deque


BASIS = (8, 4, 2, 1)
REQUIRED = (0b1011, 0b1101, 0b1110, 0b0110)
FLAG = {mask: 1 << i for i, mask in enumerate(REQUIRED)}
ALL_FLAGS = 15
EDGES = [(c, t) for c in range(4) for t in range(4) if c != t]


def step(rows, flags, edge):
    control, target = edge
    output = list(rows)
    output[target] ^= output[control]
    output = tuple(output)
    return output, flags | FLAG.get(output[target], 0)


def reconstruct(previous, state):
    edges = []
    while previous[state] is not None:
        state, edge = previous[state]
        edges.append(edge)
    return edges[::-1]


def main():
    start = (BASIS, 0)
    previous = {start: None}
    frontier = [start]
    layer_sizes = [1]
    for depth in range(6):
        next_frontier = []
        for rows, flags in frontier:
            for edge in EDGES:
                state = step(rows, flags, edge)
                if state not in previous:
                    previous[state] = ((rows, flags), edge)
                    next_frontier.append(state)
        frontier = next_frontier
        layer_sizes.append(len(frontier))
    two_cnot_endpoints = {}
    for first in EDGES:
        for second in EDGES:
            rows, flags = step(BASIS, 0, first)
            rows, flags = step(rows, flags, second)
            two_cnot_endpoints.setdefault(rows, []).append([first, second])
    found = []
    all_masks_states = 0
    for rows, flags in frontier:
        if flags != ALL_FLAGS:
            continue
        all_masks_states += 1
        if rows in two_cnot_endpoints:
            found.append({"prefix_edges": reconstruct(previous, (rows, flags)),
                          "endpoint_rows": rows,
                          "endpoint_cnot_factorizations": two_cnot_endpoints[rows]})
    print(json.dumps({"depth": 6, "layer_sizes": layer_sizes,
                      "states_at_depth_six_with_all_masks": all_masks_states,
                      "two_cnot_endpoint_maps": len(two_cnot_endpoints),
                      "found": found}, indent=2))


if __name__ == "__main__":
    main()
