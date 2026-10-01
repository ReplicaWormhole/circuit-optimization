"""Exact exploratory Schmidt-rank screen for V4 eigenspaces and six-CNOT cuts.

This script reports algebraic conditions only.  It does not infer a crossing
bound unless the symbolic result proves one.
"""

from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
EDGES = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
CUTS = ((0, 1), (0, 2), (0, 3))


def shift_matrix():
    v = sp.zeros(16)
    for state in range(16):
        v[((state & 1) << 3) | (state >> 1), state] = 1
    return v


def matrix_across_cut(vector, left):
    right = tuple(j for j in range(4) if j not in left)
    out = sp.zeros(4)
    for state in range(16):
        bits = tuple((state >> (3 - j)) & 1 for j in range(4))
        row = 2 * bits[left[0]] + bits[left[1]]
        col = 2 * bits[right[0]] + bits[right[1]]
        out[row, col] = vector[state]
    return out


def rank_two_ideal(mats):
    xs = sp.symbols('x:' + str(len(mats)))
    matrix = sum((x * m for x, m in zip(xs, mats)), sp.zeros(4))
    minors = [sp.factor(matrix.extract(rows, cols).det())
              for rows in combinations(range(4), 3)
              for cols in combinations(range(4), 3)]
    nonzero = sorted(set(p for p in minors if p != 0), key=str)
    gb = sp.groebner(nonzero, *xs, order='grevlex', domain=sp.QQ_I)
    return {
        'variables': [str(x) for x in xs],
        'basis_matrices': [[[str(z) for z in row] for row in m.tolist()] for m in mats],
        'nonzero_minor_count': len(nonzero),
        'groebner_basis': [str(p.as_expr()) for p in gb.polys],
    }


def main():
    v = shift_matrix()
    eye = sp.eye(16)
    # M_k=4P_k.  Eigenvalue i sector has dimension three.
    m1 = eye - sp.I * v - v**2 + sp.I * v**3
    assert sp.trace(m1) == 12
    columns = m1.columnspace()
    assert len(columns) == 3
    eigenspace = {}
    for left in CUTS:
        mats = [matrix_across_cut(col, left) for col in columns]
        eigenspace[''.join(map(str, left))] = rank_two_ideal(mats)
    result = {'i_eigenspace_dimension': len(columns),
              'rank_at_most_two_conditions_by_cut': eigenspace}
    path = ROOT / 'global_lower_bound7_cut_schmidt_result.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    for cut, data in eigenspace.items():
        print(cut, data['groebner_basis'])


if __name__ == '__main__':
    main()
