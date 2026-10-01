"""Exact sparse output gauge giving Schmidt rank eight on tail cut 03|12.

All identities are checked over Q(zeta_96), including gauge unitarity,
V4 diagonalization after the fixed six-CNOT prefix, and operator ranks.
The gauge is an output eigenspace operation, not a seven-CNOT synthesis.
"""

import json
from pathlib import Path

from exact_check import exact_eigenvalue_labels
from search13_crosspair_single_reroute_rank import exact_rank, realignment, CUTS
from search13_fulltail_invariant_commutant import prefix_action
from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_secondpair_mask_transport import CyclotomicPair

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_fulltail_invariant_exact_witness_result.json"


def gauge_matrix(exact):
    zero, one, imag = exact.zero, exact.one, exact.imag
    h = (exact.z**12 + exact.z**-12) * exact.half
    assert h*h + h*h == one
    g = [[zero for _ in range(16)] for _ in range(16)]
    entries = {
        0: {0: one},
        1: {1: imag*h, 14: h},
        2: {2: one},
        3: {13: one},
        4: {6: -h, 9: h},
        5: {5: h, 10: h},
        6: {4: one},
        7: {7: one},
        8: {15: one},
        9: {6: h, 9: h},
        10: {5: -h, 10: h},
        11: {11: one},
        12: {12: one},
        13: {3: one},
        14: {1: imag*h, 14: -h},
        15: {8: one},
    }
    for row, columns in entries.items():
        for col, value in columns.items():
            g[row][col] = value
    conjugates = {zero: zero, one: one, -one: -one,
                  h: h, -h: -h, imag*h: -imag*h, -imag*h: imag*h}
    gdag = [[conjugates[g[col][row]] for col in range(16)]
            for row in range(16)]
    assert exact.matmul(g, gdag) == exact.identity(16)
    return g, entries


def prefix_matrix(exact):
    perm, exponents = prefix_action()
    p = [[exact.zero for _ in range(16)] for _ in range(16)]
    for col, row in enumerate(perm):
        p[row][col] = exact.z ** (6 * exponents[col])
    return p


def main():
    source = json.loads(SOURCE.read_text())
    exact_source = exact_eigenvalue_labels(source)
    assert exact_source["exact_diagonalizer"]
    exact = CyclotomicPair()
    g, entries = gauge_matrix(exact)
    labels = exact_source["output_labels"]
    assert all(labels[row] == labels[col]
               for row, columns in entries.items() for col in columns)
    tail = tail_matrix(exact)
    gauged_tail = exact.matmul(g, tail)
    ranks = {"".join(map(str, cut)):
             exact_rank(realignment(gauged_tail, cut, exact), exact)
             for cut in CUTS}
    assert ranks["03"] == 8

    full = exact.matmul(gauged_tail, prefix_matrix(exact))
    roots = (exact.one, exact.imag, -exact.one, -exact.imag)
    shift = lambda x: ((x & 1) << 3) | (x >> 1)
    exact_diagonalizer = all(
        full[row][shift(x)] == roots[labels[row]] * full[row][x]
        for row in range(16) for x in range(16))
    assert exact_diagonalizer

    # The numerical right-kernel ansatz is now an exact finite certificate.
    re = realignment(gauged_tail, (0, 3), exact)
    zero_columns = [col for col in range(16)
                    if all(re[row][col] == exact.zero for row in range(16))]
    pair_relation = all(re[row][2] == exact.imag*re[row][1]
                        for row in range(16))
    assert zero_columns == [3, 7, 8, 9, 10, 12, 15] and pair_relation
    result = {
        "field": "Q(zeta_96)", "source": SOURCE.name,
        "prefix_cnot_count": 6, "original_tail_cnot_count": 8,
        "gauge_sparse_rows": {
            str(row): {str(col): str(exact.field.to_sympy(value))
                       for col, value in columns.items()}
            for row, columns in entries.items()},
        "gauge_unitary_exact": True,
        "preserves_output_eigenvalue_labels_exact": True,
        "diagonalizes_V4_exact": exact_diagonalizer,
        "operator_schmidt_ranks": ranks,
        "cut03_exact_zero_columns": zero_columns,
        "cut03_exact_column2_equals_i_column1": pair_relation,
        "scope": "exact rank-eight gauged eight-CNOT tail; seven-CNOT synthesis open",
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"exact_diagonalizer": exact_diagonalizer,
                      "ranks": ranks, "zero_columns": zero_columns}), flush=True)


if __name__ == "__main__":
    main()
