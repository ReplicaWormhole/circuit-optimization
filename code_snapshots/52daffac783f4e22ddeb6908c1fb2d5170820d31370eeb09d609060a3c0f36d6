"""Exact local-equivalence screen for an output CNOT on a terminal matchgate.

The final two disjoint blocks in the canonical 14-CNOT circuit are
F = exp(i*pi*(XX+YY)/8). Any computational-basis permutation of output rows
preserves diagonalization. This tests whether appending one CNOT on either
terminal pair makes that pair locally equivalent to a one-CNOT gate. If so,
the pair could be re-synthesized with one CNOT and total count would be 14.
It does not test inter-pair output permutations or more general resynthesis.
"""

import json

import sympy as s


I = s.I
R2 = s.sqrt(2)
F = s.Matrix([[1, 0, 0, 0],
              [0, 1/R2, I/R2, 0],
              [0, I/R2, 1/R2, 0],
              [0, 0, 0, 1]])
CX01 = s.Matrix([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]])
CX10 = s.Matrix([[1, 0, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0],
                 [0, 1, 0, 0]])

# Columns are an orthonormal magic basis. For a local SU(2)xSU(2) gate,
# Q^dagger U Q is real orthogonal. Eigenvalue multiplicities of M are
# invariant under arbitrary local unitaries on both sides of U.
Q = s.Matrix([[1, I, 0, 0],
              [0, 0, I, 1],
              [0, 0, I, -1],
              [1, -I, 0, 0]]) / R2


def invariant(u):
    b = Q.conjugate().T * u * Q
    m = b.T * b
    z = s.symbols("z")
    polynomial = s.factor(m.charpoly(z).as_expr(), extension=R2)
    roots = s.roots(polynomial, z)
    return str(polynomial), {str(k): v for k, v in roots.items()}


def main():
    cases = {"F": F, "CX01": CX01, "CX10": CX10,
             "CX01_F": CX01 * F, "CX10_F": CX10 * F,
             "F_CX01": F * CX01, "F_CX10": F * CX10}
    rows = {}
    for name, unitary in cases.items():
        poly, roots = invariant(unitary)
        rows[name] = {"magic_invariant_characteristic": poly,
                      "eigenvalue_multiplicities": roots,
                      "at_most_two_distinct": len(roots) <= 2}
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
