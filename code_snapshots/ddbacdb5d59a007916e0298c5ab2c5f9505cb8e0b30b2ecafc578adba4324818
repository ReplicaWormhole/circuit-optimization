"""Exact cyclotomic test of the sparse numerical rank-eight gauge ansatz.

The run254 right-kernel projector suggests seven zero realignment columns
and one exact pair relation. These are linear in each output-gauge row.
This script computes their rowwise nullspaces over Q(zeta_96), restricted
to the stable numerical support masks, without promoting a near-zero SVD
to an exact rank claim.
"""

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi
from search13_secondpair_mask_transport import CyclotomicPair

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_fulltail_invariant_exact_ansatz_result.json"
LABELS = (0, 3, 1, 1, 2, 0, 2, 2, 0, 2, 0, 0, 3, 1, 3, 0)
ZERO_COLUMNS = (3, 7, 8, 9, 10, 12, 15)
PAIR_RELATION = (2, 1, "i")  # column 2 minus i times column 1
MASKS = (
    (0, 5, 10, 15), (1, 14), (2, 13), (2, 13), (6, 9),
    (5, 10), (4, 7), (4, 7), (0, 5, 10, 15), (6, 9),
    (0, 5, 8, 10, 11, 15), (0, 5, 8, 10, 11, 15),
    (12,), (3,), (1, 14), (8, 11),
)


def tail_matrix(exact):
    gates = json.loads(SOURCE.read_text())["gates"][13:]
    out = exact.identity(16)
    cnot_count = 0
    for gate in gates:
        name = gate["gate"]
        if name == "cx":
            cnot_count += 1
            cmask, tmask = 1 << (3-gate["control"]), 1 << (3-gate["target"])
            next_rows = [None] * 16
            for y in range(16):
                z = y ^ tmask if y & cmask else y
                next_rows[z] = out[y]
            out = next_rows
        else:
            theta = rational_pi(gate["theta"])
            if name == "rz":
                local = exact.rz(theta)
            elif name in ("rx", "ry"):
                c, s = exact.trig(theta / 2)
                if name == "rx":
                    local = ((c, -exact.imag*s), (-exact.imag*s, c))
                else:
                    local = ((c, -s), (s, c))
            else:
                raise ValueError(f"unexpected tail gate {name}")
            bit = 1 << (3-gate["qubit"])
            for y in range(16):
                if y & bit:
                    continue
                z = y ^ bit
                a, b = out[y], out[z]
                out[y] = [local[0][0]*u + local[0][1]*v
                          for u, v in zip(a, b)]
                out[z] = [local[1][0]*u + local[1][1]*v
                          for u, v in zip(a, b)]
    assert cnot_count == 8
    return out


def right_bits(index):
    return (((index >> 2) & 1) << 1) | ((index >> 1) & 1)


def left_bits(index):
    return (((index >> 3) & 1) << 1) | (index & 1)


def nullspace(rows, ncols, exact):
    a = [row[:] for row in rows]
    pivots = []
    for col in range(ncols):
        pivot = next((i for i in range(len(pivots), len(a))
                      if a[i][col] != exact.zero), None)
        if pivot is None:
            continue
        pos = len(pivots)
        a[pos], a[pivot] = a[pivot], a[pos]
        scalar = a[pos][col]
        a[pos] = [v / scalar for v in a[pos]]
        for i in range(len(a)):
            if i == pos or a[i][col] == exact.zero:
                continue
            scalar = a[i][col]
            a[i] = [v - scalar*w for v, w in zip(a[i], a[pos])]
        pivots.append(col)
    free = [j for j in range(ncols) if j not in pivots]
    basis = []
    for col in free:
        vector = [exact.zero] * ncols
        vector[col] = exact.one
        for row, pivot in enumerate(pivots):
            vector[pivot] = -a[row][col]
        basis.append(vector)
    return basis


def constraints(tail, y, mask, exact):
    yr = right_bits(y)
    rows = []
    for c in ZERO_COLUMNS:
        if c // 4 != yr:
            continue
        xr = c % 4
        for x in range(16):
            if right_bits(x) == xr:
                rows.append([tail[s][x] for s in mask])
    if yr == 0:
        for xl in range(4):
            x2 = next(x for x in range(16) if right_bits(x) == 2
                      and left_bits(x) == xl)
            x1 = next(x for x in range(16) if right_bits(x) == 1
                      and left_bits(x) == xl)
            rows.append([tail[s][x2] - exact.imag*tail[s][x1] for s in mask])
    return rows


def main():
    exact = CyclotomicPair()
    tail = tail_matrix(exact)
    sectors = {k: [i for i, label in enumerate(LABELS) if label == k]
               for k in range(4)}
    results = []
    for y, mask in enumerate(MASKS):
        assert all(s in sectors[LABELS[y]] for s in mask)
        equations = constraints(tail, y, mask, exact)
        basis = nullspace(equations, len(mask), exact)
        results.append({"row": y, "eigenvalue_power": LABELS[y],
                        "allowed_columns": mask, "equation_count": len(equations),
                        "exact_nullity": len(basis),
                        "basis": [[str(exact.field.to_sympy(z)) for z in v]
                                  for v in basis]})
        print(json.dumps({k: results[-1][k] for k in ("row", "allowed_columns",
              "exact_nullity")}), flush=True)
    result = {"field": "Q(zeta_96)", "source": SOURCE.name,
              "zero_realign_columns": ZERO_COLUMNS,
              "pair_relation": PAIR_RELATION, "rows": results,
              "scope": "exact linear ansatz feasibility only; row unitarity not yet solved"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
