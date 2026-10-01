"""Exact cut-rank screen for one final-layer CNOT rerouted to 03 or 12.

For 196 five-prefix near cases, hold the exact first 02/13 operator and the
resulting target four-CNOT suffix fixed. Replace exactly one CNOT in the
chronological 01,01,23,23 suffix by a directed 03 or 12 CNOT. Arbitrary
one-qubit gates are allowed around all four CNOTs. Across any cut, k crossing
CNOTs give operator Schmidt rank at most 2**k. Compute target ranks exactly
over Q(zeta_96), checking single-wire cuts before balanced cuts.
"""

import json
from collections import Counter
from fractions import Fraction

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair, PAIRS, enumerate_near_cases
from search13_crosspair_second_layer_rank import CENTER, four_wire_factor


CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
BASE_EDGES = ((0, 1), (0, 1), (2, 3), (2, 3))


def realignment(operator, cut, exact):
    left = tuple(cut)
    right = tuple(wire for wire in range(4) if wire not in left)
    dl, dr = 1 << len(left), 1 << len(right)
    out = [[exact.zero for _ in range(dr*dr)] for _ in range(dl*dl)]

    def coordinate(basis, wires):
        value = 0
        for wire in wires:
            value = (value << 1) | ((basis >> (3-wire)) & 1)
        return value

    for y in range(16):
        yl, yr = coordinate(y, left), coordinate(y, right)
        for x in range(16):
            xl, xr = coordinate(x, left), coordinate(x, right)
            out[yl*dl+xl][yr*dr+xr] = operator[y][x]
    return out


def exact_rank(matrix, exact):
    a = [row[:] for row in matrix]
    nrows, ncols = len(a), len(a[0])
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, nrows) if a[r][col] != exact.zero), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        denominator = a[rank][col]
        for row in range(rank+1, nrows):
            if a[row][col] == exact.zero:
                continue
            factor = a[row][col] / denominator
            for j in range(col, ncols):
                a[row][j] -= factor * a[rank][j]
        rank += 1
        if rank == nrows:
            break
    return rank


def topologies():
    rows = []
    for position in range(4):
        for new_pair in ((0, 3), (1, 2)):
            for direction in (new_pair, new_pair[::-1]):
                edges = list(BASE_EDGES)
                edges[position] = direction
                crossing = {"".join(map(str, cut)): sum(
                    (a in cut) != (b in cut) for a, b in edges)
                    for cut in CUTS}
                rows.append({"replaced_position": position,
                             "new_directed_edge": list(direction),
                             "edges": [list(edge) for edge in edges],
                             "crossing_counts": crossing})
    assert len(rows) == 16
    return rows


def main():
    cases = enumerate_near_cases()
    assert len(cases) == 196
    unique = {}
    for case in cases:
        key = (tuple(case["endpoint"]), case["missing_mask"],
               int(case["missing_coefficient_pi_over_8"]) % 16)
        unique.setdefault(key, 0)
        unique[key] += 1
    assert len(unique) == 51

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
    baseline_second = exact.matmul(exact.kron(f2, f2), center)

    routes = topologies()
    phase_cache = {}
    rank_histogram = {"".join(map(str, cut)): Counter() for cut in CUTS}
    result_rows = []
    survivors = []
    for (endpoint, mask, coefficient), multiplicity in unique.items():
        key = endpoint, mask
        if key not in phase_cache:
            bits = exact.mask_bits(endpoint, mask)
            q02 = exact.transported_factor(endpoint, PAIRS[0], bits)
            q13 = exact.transported_factor(endpoint, PAIRS[1], bits)
            q16 = four_wire_factor(q02, q13, exact)
            phase_cache[key] = bits, exact.matmul(baseline_second, q16)
        bits, bq = phase_cache[key]
        angle_cos, angle_sin = exact.trig(Fraction(coefficient, 16))
        target = exact.add_scaled(angle_cos, baseline_second,
                                  -imag*angle_sin, bq)
        ranks = {}
        for cut in CUTS[:4]:
            label = "".join(map(str, cut))
            ranks[label] = exact_rank(realignment(target, cut, exact), exact)
            rank_histogram[label][ranks[label]] += 1
        single_survivors = [route for route in routes if all(
            ranks[str(wire)] <= 2**route["crossing_counts"][str(wire)]
            for wire in range(4))]
        if single_survivors:
            for cut in CUTS[4:]:
                label = "".join(map(str, cut))
                ranks[label] = exact_rank(realignment(target, cut, exact), exact)
                rank_histogram[label][ranks[label]] += 1
        compatible = [route for route in single_survivors if all(
            ranks["".join(map(str, cut))] <=
            2**route["crossing_counts"]["".join(map(str, cut))]
            for cut in CUTS[4:])]
        row = {"endpoint": list(endpoint), "missing_mask": mask,
               "coefficient_pi_over_8_mod16": coefficient,
               "post_prefix_wire_mask": bits,
               "multiplicity": multiplicity,
               "target_schmidt_ranks": ranks,
               "single_cut_surviving_route_count": len(single_survivors),
               "all_cut_surviving_route_count": len(compatible)}
        result_rows.append(row)
        for route in compatible:
            survivors.append({**row, "route": route})

    result = {"near_cases": len(cases), "distinct_targets": len(unique),
              "distinct_transported_generators": len(phase_cache),
              "topology_count": len(routes),
              "topologies": routes,
              "rank_histogram": {key: dict(value) for key, value in rank_histogram.items()},
              "compatible_case_route_count": len(survivors),
              "compatible_case_routes": survivors,
              "target_rows": result_rows}
    out = ROOT / "search13_crosspair_single_reroute_rank_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("topologies", "compatible_case_routes", "target_rows")},
                     indent=2))


if __name__ == "__main__":
    main()
