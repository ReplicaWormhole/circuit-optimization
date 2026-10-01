"""Numerical balanced-cut rank screen for two complex eigenspace rotations.

The fixed exact14 prefix is six CNOTs.  Gauge the full eight-CNOT tail and
the fixed final four-CNOT layer by two disjoint, same-eigenvalue row SU(2)
rotations.  Rank-compatible three- or seven-CNOT graphs are necessary only.
"""

import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
import json

import numpy as np

from delete14_search import fixed_matrix
from search13_joint_gauge_fit import BASE, ROOT
from search13_output_eigenrow_gauge_rank import LABELS


CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
EDGES = tuple(combinations(range(4), 2))
PHASE_NUMERATORS = (-2, -1, 1, 2)  # multiples of pi/4, all non-real
ANGLE_NUMERATORS = (1, 2)           # multiples of pi/8


def target(gates):
    out = np.eye(16, dtype=complex)
    for gate in gates:
        out = fixed_matrix(gate) @ out
    return out


def realigned_rank(matrix, cut):
    tensor = matrix.reshape((2,)*8)
    left = tuple(cut) + tuple(4+j for j in cut)
    right = tuple(j for j in range(8) if j not in left)
    dim = 1 << (2*len(cut))
    s = np.linalg.svd(tensor.transpose(left+right).reshape(dim, 256//dim),
                      compute_uv=False)
    return int(np.count_nonzero(s > 1e-9))


def ranks(matrix):
    return tuple(realigned_rank(matrix, cut) for cut in CUTS)


def graph_catalog(n):
    graphs = []
    for edge_ids in combinations_with_replacement(range(6), n):
        edges = tuple(EDGES[i] for i in edge_ids)
        caps = tuple(2**sum((a in cut) != (b in cut) for a, b in edges)
                     for cut in CUTS)
        graphs.append((edges, caps))
    return graphs


def gauged(matrix, pair1, pair2, angle1, angle2, phase1, phase2):
    out = matrix.copy()
    for (first, second), angle, phase in ((pair1, angle1, phase1),
                                         (pair2, angle2, phase2)):
        a, b = out[first].copy(), out[second].copy()
        co, si = np.cos(angle), np.sin(angle)
        ep = np.exp(1j*phase)
        out[first] = co*a - np.conj(ep)*si*b
        out[second] = ep*si*a + co*b
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()
    gates = json.loads(BASE.read_text())["gates"]
    full = target(gates[13:])
    final = target(gates[39:])
    assert ranks(full)[4:] == (16, 16, 16)
    assert ranks(final)[4:] == (1, 16, 16)
    pairs = [pair for pair in combinations(range(16), 2)
             if LABELS[pair[0]] == LABELS[pair[1]]]
    graphs3, graphs7 = graph_catalog(3), graph_catalog(7)
    rows = []
    hist_final = Counter()
    hist_full = Counter()
    best_final = []
    best_full = []
    total = 0
    for first, second in combinations(pairs, 2):
        if set(first) & set(second):
            continue
        for angles in product(ANGLE_NUMERATORS, repeat=2):
            for phases in product(PHASE_NUMERATORS, repeat=2):
                total += 1
                argspec = (first, second, *(a*np.pi/8 for a in angles),
                           *(p*np.pi/4 for p in phases))
                fr = ranks(gauged(final, *argspec))
                tr = ranks(gauged(full, *argspec))
                hist_final[str(fr)] += 1
                hist_full[str(tr)] += 1
                c3 = [i for i, (_, cap) in enumerate(graphs3)
                      if all(r <= c for r, c in zip(fr, cap))]
                c7 = [i for i, (_, cap) in enumerate(graphs7)
                      if all(r <= c for r, c in zip(tr, cap))]
                row = {"pairs": [list(first), list(second)],
                       "angle_pi_over_8": angles,
                       "phase_pi_over_4": phases,
                       "final_ranks": fr, "full_ranks": tr,
                       "compatible_three_graphs": c3,
                       "compatible_seven_graph_count": len(c7),
                       "best_seven_graph": c7[0] if c7 else None}
                if c3:
                    rows.append(row)
                best_final.append(row)
                best_full.append(row)
    best_final.sort(key=lambda r: (sum(r["final_ranks"][4:]),
                                   r["final_ranks"][4:], r["pairs"]))
    best_full.sort(key=lambda r: (sum(r["full_ranks"][4:]),
                                  r["full_ranks"][4:], r["pairs"]))
    result = {"base": BASE.name, "total_gauges": total,
              "disjoint_eigenrow_pairs": sum(not(set(a) & set(b))
                  for a, b in combinations(pairs, 2)),
              "angles_pi_over_8": ANGLE_NUMERATORS,
              "phases_pi_over_4": PHASE_NUMERATORS,
              "rank_tolerance": 1e-9,
              "three_cnot_graphs": len(graphs3),
              "seven_cnot_graphs": len(graphs7),
              "three_compatible_count": len(rows),
              "three_compatible": rows,
              "top_final": best_final[:args.top],
              "top_full": best_full[:args.top],
              "final_rank_histogram": dict(hist_final),
              "full_rank_histogram": dict(hist_full),
              "graph3": [[list(edge) for edge in graph[0]] for graph in graphs3],
              "graph7": [[list(edge) for edge in graph[0]] for graph in graphs7]}
    path = ROOT / "search13_joint_gauge_phased_rank_result.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"result": path.name, "total_gauges": total,
                      "three_compatible_count": len(rows),
                      "best_final": best_final[0],
                      "best_full": best_full[0]}, sort_keys=True))


if __name__ == "__main__":
    main()
