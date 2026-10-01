"""Direction-aware CNOT tensor mincuts for eight near-tail schedules.

Each CNOT is factored exactly as sum_c |c><c| tensor X^c, a rank-two
bond between a control factor and a target factor. All wire segments
also have dimension two. Arbitrary one-qubit gates sit on wire segments.
For every orientation assignment, graph mincuts bound realignment ranks
of the exact rank-eight-gauged target. Only masks 83 and 163 are used here;
their target ranks are certified modulo 97 and 193 in run 284.
"""

from itertools import product
import json
from pathlib import Path

import networkx as nx


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_rank8_neighborhood_mixedcut_result.json"
OUT = ROOT / "search13_rank8_directional_mincut_result.json"
MASKS = (83, 163)
INF = 100


def graph(schedule, directions):
    graph = nx.DiGraph()

    def bond(a, b):
        graph.add_edge(a, b, capacity=1)
        graph.add_edge(b, a, capacity=1)

    factors = {}
    for index, ((a, b), direction) in enumerate(zip(schedule, directions)):
        c, t = (a, b) if direction == 0 else (b, a)
        control_node, target_node = 8 + 2*index, 9 + 2*index
        factors[(index, c)] = control_node
        factors[(index, t)] = target_node
        bond(control_node, target_node)
    for wire in range(4):
        visits = [factors[(index, wire)] for index, pair in enumerate(schedule)
                  if wire in pair]
        path = [4 + wire] + visits + [wire]
        for a, b in zip(path, path[1:]):
            bond(a, b)
    return graph


def mincut(g, mask):
    augmented = g.copy()
    source, sink = 30, 31
    for boundary in range(8):
        if mask & (1 << boundary):
            augmented.add_edge(source, boundary, capacity=INF)
        else:
            augmented.add_edge(boundary, sink, capacity=INF)
    value, _ = nx.minimum_cut(augmented, source, sink,
                              flow_func=nx.algorithms.flow.preflow_push)
    return int(value)


def main():
    data = json.loads(SOURCE.read_text())
    ranks = {mask: data["target_modular_rank_lower_bounds"][str(mask)]
             for mask in MASKS}
    assert ranks == {83: 14, 163: 14}
    rows = []
    for schedule in data["surviving_schedules"]:
        acceptable = []
        minima = {mask: 99 for mask in MASKS}
        maxima = {mask: 0 for mask in MASKS}
        for directions in product(range(2), repeat=7):
            g = graph(schedule, directions)
            cuts = {mask: mincut(g, mask) for mask in MASKS}
            for mask in MASKS:
                minima[mask] = min(minima[mask], cuts[mask])
                maxima[mask] = max(maxima[mask], cuts[mask])
            if all(ranks[mask] <= 2**cuts[mask] for mask in MASKS):
                acceptable.append({"directions": list(directions),
                                   "mincuts": {str(k): v for k, v in cuts.items()}})
        rows.append({"schedule": schedule,
                     "tested_orientations": 128,
                     "orientations_surviving_two_cuts": len(acceptable),
                     "example_survivor": acceptable[0] if acceptable else None,
                     "minimum_mincuts": {str(k): v for k, v in minima.items()},
                     "maximum_mincuts": {str(k): v for k, v in maxima.items()}})
    result = {"source": SOURCE.name, "target_rank_lower_bounds": ranks,
              "tested_boundary_masks": MASKS,
              "rows": rows,
              "schedule_exclusions": sum(row["orientations_surviving_two_cuts"] == 0
                                         for row in rows),
              "scope": "exact necessary directional tensor-network test for eight listed schedules and two cuts only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"schedule_exclusions": result["schedule_exclusions"],
                      "surviving_orientations_per_schedule": [
                          row["orientations_surviving_two_cuts"] for row in rows]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
