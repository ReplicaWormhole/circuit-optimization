"""Exact parity-conditioned fermionic Fourier diagonalizer of qubit V4.

In N-particle sector, the qubit cyclic shift is fermionic mode translation
with boundary condition eta=(-1)^(N-1): periodic for odd N, antiperiodic
for even N. Use the corresponding one-particle Fourier matrix F_N and
its exterior powers. The resulting 16x16 Fock transform is checked
entrywise over Q(zeta_8) for U V4 = D U.
"""

import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp

from algebraic14_orbit_shear import rotate


ROOT = Path(__file__).parent
FIELD = sp.QQ.cyclotomic_field(8)
ZETA = FIELD.unit
ZERO, ONE = FIELD.zero, FIELD.one
HALF = FIELD.from_sympy(sp.Rational(1, 2))
ROOTS = (ONE, ZETA ** 2, -ONE, -(ZETA ** 2))


def occupied(x):
    return tuple(j for j in range(4) if x & (1 << (3 - j)))


def determinant(matrix):
    n = len(matrix)
    total = ZERO
    for permutation in itertools.permutations(range(n)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(n) for j in range(i + 1, n))
        term = ONE
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        total += -term if inversions & 1 else term
    return total


def single_particle(n):
    offset = 0 if n & 1 else 1
    return [[HALF * ZETA ** ((2 * k + offset) * j)
             for j in range(4)] for k in range(4)]


def fock_transform():
    unitary = [[ZERO for _ in range(16)] for _ in range(16)]
    for output in range(16):
        out_modes = occupied(output)
        n = len(out_modes)
        f = single_particle(n)
        for source in range(16):
            in_modes = occupied(source)
            if len(in_modes) != n:
                continue
            minor = [[f[k][j] for j in in_modes] for k in out_modes]
            unitary[output][source] = determinant(minor)
    return unitary


def main():
    u = fock_transform()
    labels = []
    for row in u:
        matching = [index for index, root in enumerate(ROOTS)
                    if all(row[rotate(x)] == root * row[x]
                           for x in range(16))]
        assert len(matching) == 1, matching
        labels.append(matching[0])
    numeric = np.array([[complex(FIELD.to_sympy(value).evalf())
                         for value in row] for row in u])
    unitary_error = float(np.max(np.abs(
        numeric @ numeric.conj().T - np.eye(16))))
    result = {"cyclotomic_conductor": 8,
              "exact_cycle_diagonalizer": True,
              "output_labels": labels,
              "multiplicities": [labels.count(k) for k in range(4)],
              "numeric_unitarity_error": unitary_error,
              "boundary_condition_by_particle_number": {
                  str(n): "periodic" if n & 1 else "antiperiodic"
                  for n in range(5)}}
    (ROOT / "algebraic14_parity_fermion_fourier_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
