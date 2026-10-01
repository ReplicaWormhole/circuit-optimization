"""Exact restricted obstruction for the closest seven-CNOT tail skeleton.

Across the mixed boundary cut {out0,out1,in1,in3}, the ordered skeleton
02,02,03,13,01,01,23 has a two-bond tensor-network cut and hence rank
at most four for arbitrary local gates and either CNOT direction. A 5x5
minor of the exact rank-eight-gauged target is nonzero in Q(zeta_96).
"""

import json
from pathlib import Path

import sympy as sp

from search13_crosspair_single_reroute_rank import exact_rank
from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_fulltail_invariant_exact_witness import gauge_matrix
from search13_fulltail_rank8_fit import SCHEDULES
from search13_secondpair_mask_transport import CyclotomicPair

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_fulltail_rank8_mixedcut_certificate_result.json"
SCHEDULE = SCHEDULES[0]
MASK = 163
LEFT = (0, 1, 5, 7)
RIGHT = tuple(j for j in range(8) if j not in LEFT)
MINOR_ROWS = (9, 1, 6, 7, 15)
MINOR_COLS = (8, 14, 5, 15, 13)


def boundary_coordinate(bits, axes):
    value = 0
    for axis in axes:
        value = 2 * value + bits[axis]
    return value


def realignment(matrix, exact):
    out = [[exact.zero for _ in range(16)] for _ in range(16)]
    for y in range(16):
        for x in range(16):
            bits = [(y >> (3-j)) & 1 for j in range(4)] + [
                    (x >> (3-j)) & 1 for j in range(4)]
            out[boundary_coordinate(bits, LEFT)][
                boundary_coordinate(bits, RIGHT)] = matrix[y][x]
    return out


def determinant(matrix, exact):
    a = [row[:] for row in matrix]
    det = exact.one
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a))
                      if a[i][col] != exact.zero), None)
        if pivot is None:
            return exact.zero
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        diagonal = a[col][col]
        det *= diagonal
        for row in range(col + 1, len(a)):
            if a[row][col] == exact.zero:
                continue
            factor = a[row][col] / diagonal
            for j in range(col, len(a)):
                a[row][j] -= factor * a[col][j]
    return det


def wire_graph():
    edges = []
    for wire in range(4):
        gates = [8+t for t, pair in enumerate(SCHEDULE) if wire in pair]
        path = [4+wire] + gates + [wire]
        edges.extend((a, b, wire) for a, b in zip(path, path[1:]))
    return edges


def graph_cut():
    edges = wire_graph()
    best = None
    for assignments in range(1 << 7):
        sides = [(MASK >> j) & 1 for j in range(8)]
        sides += [(assignments >> j) & 1 for j in range(7)]
        crossing = [(a, b, wire) for a, b, wire in edges
                    if sides[a] != sides[b]]
        if best is None or len(crossing) < len(best[1]):
            best = (assignments, crossing)
    assert len(best[1]) == 2
    return best


def main():
    exact = CyclotomicPair()
    gauge, _ = gauge_matrix(exact)
    target = exact.matmul(gauge, tail_matrix(exact))
    matrix = realignment(target, exact)
    sub = [[matrix[row][col] for col in MINOR_COLS]
           for row in MINOR_ROWS]
    det = determinant(sub, exact)
    expected = exact.field.from_sympy(sp.Rational(3, 16)) * exact.z**32
    assert det == expected and det != exact.zero
    rank = exact_rank(matrix, exact)
    assert rank == 14
    assignment, crossing = graph_cut()
    result = {
        "field": "Q(zeta_96)",
        "fixed_prefix": "exact14 first six CNOTs",
        "target": "exact rank-eight output gauge times exact14 eight-CNOT tail",
        "ordered_unoriented_tail_skeleton": [list(edge) for edge in SCHEDULE],
        "boundary_axis_order": "out0,out1,out2,out3,in0,in1,in2,in3",
        "left_boundary_axes": list(LEFT),
        "right_boundary_axes": list(RIGHT),
        "boundary_cut_mask": MASK,
        "gate_vertex_side_mask": assignment,
        "cut_wire_segment_edges": [{"first_vertex": a,
                                    "second_vertex": b, "wire": wire}
                                   for a, b, wire in crossing],
        "tensor_network_mincut_bonds": 2,
        "skeleton_rank_upper_bound": 4,
        "target_exact_rank": rank,
        "nonzero_minor_rows": list(MINOR_ROWS),
        "nonzero_minor_columns": list(MINOR_COLS),
        "nonzero_5x5_minor_determinant": "3*zeta_96^32/16",
        "exactly_excludes_this_ordered_unoriented_skeleton": True,
        "scope": "arbitrary one-qubit gates and both directions for these seven physical pairs in this order only; other tail schedules untested",
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"target_rank": rank, "rank_cap": 4,
                      "minor": result["nonzero_5x5_minor_determinant"],
                      "cut_edges": result["cut_wire_segment_edges"]}), flush=True)


if __name__ == "__main__":
    main()
