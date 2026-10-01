"""Global seven-CNOT pair-schedule screen for one exact gauged tail.

Enumerate all 6^7 chronological unoriented pair schedules. After exact
ordinary-cut and residual-commutant necessary screens, apply the CNOT
rank-two-factor tensor mincut at mixed boundary masks 83 and 163. The
exact target has rank at least 14 at both masks, so each cut needs four
dimension-two bonds. CNOT directions do not change this abstract graph.
"""

from collections import Counter, deque
from itertools import combinations, product
import json
from pathlib import Path

from search13_fulltail_rank8_topology import CUTS, crossing, reverse_cone


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_rank8_global_twocut_result.json"
EDGES = tuple(combinations(range(4), 2))
MASKS = (83, 163)


def bond_graph(schedule):
    edges = [(8 + 2*t, 9 + 2*t) for t in range(7)]
    for wire in range(4):
        visits = [8 + 2*t + (wire == b)
                  for t, (a, b) in enumerate(schedule) if wire == a or wire == b]
        path = [4 + wire] + visits + [wire]
        edges.extend(zip(path, path[1:]))
    assert len(edges) == 25
    return edges


def edge_connectivity_up_to_four(edges, mask):
    """Unit-capacity undirected-bond maxflow, stopped at target value four."""
    residual = [dict() for _ in range(24)]
    for a, b in edges:
        assert b not in residual[a]
        residual[a][b] = 1
        residual[b][a] = 1
    source, sink = 22, 23
    for boundary in range(8):
        if mask & (1 << boundary):
            residual[source][boundary] = 4
            residual[boundary][source] = 0
        else:
            residual[boundary][sink] = 4
            residual[sink][boundary] = 0
    flow = 0
    while flow < 4:
        previous = {source: -1}
        queue = deque([source])
        while queue and sink not in previous:
            node = queue.popleft()
            for neighbor, capacity in residual[node].items():
                if capacity and neighbor not in previous:
                    previous[neighbor] = node
                    queue.append(neighbor)
        if sink not in previous:
            break
        node = sink
        while node != source:
            parent = previous[node]
            residual[parent][node] -= 1
            residual[node][parent] += 1
            node = parent
        flow += 1
    return flow


def main():
    target = json.loads((ROOT / "search13_fulltail_invariant_exact_witness_result.json").read_text())
    required = tuple(target["operator_schmidt_ranks"]["".join(map(str, cut))]
                     for cut in CUTS)
    lower = json.loads((ROOT / "search13_rank8_neighborhood_mixedcut_result.json").read_text())
    assert all(lower["target_modular_rank_lower_bounds"][str(mask)] >= 14
               for mask in MASKS)
    commutant = json.loads((ROOT / "search13_fulltail_invariant_commutant_result.json").read_text())
    forbidden = {tuple(row) for row in commutant["certified_no_commutant_subsets"]}
    counts = Counter()
    mincut_profiles = Counter()
    survivor_examples = []
    for schedule in product(EDGES, repeat=7):
        crosses = tuple(crossing(schedule, cut) for cut in CUTS)
        if any(rank > 2**count for rank, count in zip(required, crosses)):
            counts["ordinary_cut_rejected"] += 1
            continue
        if any(reverse_cone(schedule, j) in forbidden for j in range(4)):
            counts["commutant_rejected"] += 1
            continue
        graph = bond_graph(schedule)
        first = edge_connectivity_up_to_four(graph, MASKS[0])
        if first < 4:
            counts["mask83_rejected"] += 1
            mincut_profiles[(first, -1)] += 1
            continue
        second = edge_connectivity_up_to_four(graph, MASKS[1])
        mincut_profiles[(first, second)] += 1
        if second < 4:
            counts["mask163_rejected"] += 1
            continue
        counts["survives_all"] += 1
        if len(survivor_examples) < 100:
            survivor_examples.append([list(edge) for edge in schedule])
    assert sum(counts.values()) == 6**7
    result = {"target": "exact rank-eight gauged tail after fixed six-CNOT prefix",
              "target_modular_rank_lower_bound_at_masks": {"83": 14, "163": 14},
              "all_ordered_unoriented_seven_pair_schedules": 6**7,
              "exclusive_counts": dict(counts),
              "mincut_profile_counts": {f"{a},{b}": n for (a, b), n in sorted(mincut_profiles.items())},
              "up_to_100_survivor_examples": survivor_examples,
              "scope": "necessary screens for this fixed exact output gauge and six-CNOT prefix; no global 13-CNOT bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["exclusive_counts"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
