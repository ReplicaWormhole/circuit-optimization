"""Finite exact output-eigenspace gauge screen of the 14-CNOT tail.

At the exact six-CNOT prefix, the full tail has eight CNOTs; after the first
02/13 matchgate layer, the residual final layer has four. Apply a single
same-eigenvalue row swap or Givens rotation with rational-pi angle to both
fixed target operators. Exact Q(zeta_96) Schmidt ranks are screened against
all undirected seven- and three-CNOT interaction multigraphs respectively.
This is a finite necessary test, not an optimality proof or synthesis.
"""

import json
from fractions import Fraction
from itertools import combinations, combinations_with_replacement

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair
from search13_crosspair_second_layer_rank import CENTER, four_wire_factor
from search13_crosspair_single_reroute_rank import CUTS, exact_rank, realignment


LABELS = (0, 3, 1, 1, 2, 0, 2, 2, 0, 2, 0, 0, 3, 1, 3, 0)
EDGES = tuple(combinations(range(4), 2))
ANGLES = (Fraction(1, 8), Fraction(1, 4), Fraction(1, 2))


def graph_catalog(cnot_count):
    graphs = []
    for edge_ids in combinations_with_replacement(range(len(EDGES)), cnot_count):
        edges = tuple(EDGES[index] for index in edge_ids)
        capacities = {"".join(map(str, cut)): 2**sum(
            (a in cut) != (b in cut) for a, b in edges) for cut in CUTS}
        graphs.append({"edges": [list(edge) for edge in edges],
                       "capacities": capacities})
    return graphs


def gauged_rows(matrix, first, second, kind, angle, exact):
    out = [row[:] for row in matrix]
    a, b = matrix[first], matrix[second]
    if kind == "swap":
        out[first], out[second] = b[:], a[:]
    else:
        co, si = exact.trig(angle)
        out[first] = [co*x - si*y for x, y in zip(a, b)]
        out[second] = [si*x + co*y for x, y in zip(a, b)]
    return out


def ranks(matrix, exact):
    return {"".join(map(str, cut)):
            exact_rank(realignment(matrix, cut, exact), exact) for cut in CUTS}


def compatible_graphs(rank_map, graphs):
    return [index for index, graph in enumerate(graphs)
            if all(rank <= graph["capacities"][cut]
                   for cut, rank in rank_map.items())]


def main():
    exact = CyclotomicPair()
    zero, one, imag = exact.zero, exact.one, exact.imag
    x = [[zero, one], [one, zero]]
    y = [[zero, -imag], [imag, zero]]
    xx, yy = exact.kron(x, x), exact.kron(y, y)
    co, si = exact.trig(Fraction(1, 8))
    f2 = exact.matmul(exact.add_scaled(co, exact.eye4, imag*si, xx),
                      exact.add_scaled(co, exact.eye4, imag*si, yy))
    center = exact.kron(exact.kron(exact.zyz(CENTER[0]), exact.zyz(CENTER[1])),
                        exact.kron(exact.zyz(CENTER[2]), exact.zyz(CENTER[3])))
    final = exact.matmul(exact.kron(f2, f2), center)
    first = four_wire_factor(exact.f, exact.f, exact)
    input_local = exact.kron(
        exact.kron(exact.local[0], exact.local[1]),
        exact.kron(exact.local[2], exact.local[3]))
    full = exact.matmul(exact.matmul(final, first), input_local)
    graph3, graph7 = graph_catalog(3), graph_catalog(7)
    assert len(graph3) == 56 and len(graph7) == 792

    pairs = [pair for pair in combinations(range(16), 2)
             if LABELS[pair[0]] == LABELS[pair[1]]]
    assert len(pairs) == 27
    gauge_specs = [(a, b, "swap", None) for a, b in pairs]
    gauge_specs += [(a, b, "givens", angle)
                    for a, b in pairs for angle in ANGLES]
    assert len(gauge_specs) == 108

    baseline = {"final_ranks": ranks(final, exact),
                "full_ranks": ranks(full, exact)}
    rows = []
    final_compatible = []
    full_compatible = []
    for first_row, second_row, kind, angle in gauge_specs:
        gf = gauged_rows(final, first_row, second_row, kind, angle, exact)
        gt = gauged_rows(full, first_row, second_row, kind, angle, exact)
        rf, rt = ranks(gf, exact), ranks(gt, exact)
        cf, ct = compatible_graphs(rf, graph3), compatible_graphs(rt, graph7)
        record = {"rows": [first_row, second_row],
                  "eigenvalue_power": LABELS[first_row],
                  "kind": kind,
                  "angle_pi_multiple": str(angle) if angle is not None else None,
                  "final_ranks": rf, "full_ranks": rt,
                  "three_cnot_final_graph_count": len(cf),
                  "seven_cnot_full_graph_count": len(ct),
                  "three_cnot_graph_indices": cf,
                  "seven_cnot_graph_indices": ct}
        rows.append(record)
        if cf:
            final_compatible.append(record)
        if ct:
            full_compatible.append(record)

    result = {"output_eigenvalue_powers": list(LABELS),
              "within_eigenspace_pairs": len(pairs),
              "gauges": len(gauge_specs),
              "angle_pi_multiples": [str(angle) for angle in ANGLES],
              "three_cnot_graphs": len(graph3),
              "seven_cnot_graphs": len(graph7),
              "baseline": baseline,
              "three_cnot_final_compatible_count": len(final_compatible),
              "seven_cnot_full_compatible_count": len(full_compatible),
              "final_compatible_gauges": final_compatible,
              "full_compatible_gauges": full_compatible,
              "graph3": graph3, "graph7": graph7,
              "records": rows}
    out = ROOT / "search13_output_eigenrow_gauge_rank_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("graph3", "graph7", "records",
                                     "final_compatible_gauges", "full_compatible_gauges")},
                     indent=2))


if __name__ == "__main__":
    main()
