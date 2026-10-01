"""Exact finite checks supporting the global four-CNOT lower bound.

The proof is in global_lower_bound4_spectral_path.md.  This script checks the
weight-one Fourier matrix elements and enumerates four-edge cycle schedules
under the support-growth necessary condition; it does not prove five CNOTs.
"""

from fractions import Fraction
from itertools import permutations


def times_i(k):
    return ((1, 0), (0, 1), (-1, 0), (0, -1))[k % 4]


def fourier_z_entry(a, b):
    # phi_k = (1/2) sum_j i^(-kj) |e_j>; Z_0 = I - 2|e_0><e_0|.
    real = 0
    imag = 0
    for j in range(4):
        zr, zi = times_i((a - b) * j)
        sign = -1 if j == 0 else 1
        real += sign * zr
        imag += sign * zi
    return Fraction(real, 4), Fraction(imag, 4)


def backward_support(schedule, output_wire):
    support = {output_wire}
    for a, b in reversed(schedule):
        if a in support or b in support:
            support.update((a, b))
    return support


def main():
    entries = [[fourier_z_entry(a, b) for b in range(4)] for a in range(4)]
    assert all(entries[a][b] == ((Fraction(1, 2), Fraction(0))
                                  if a == b else (Fraction(-1, 2), Fraction(0)))
               for a in range(4) for b in range(4))
    assert entries[0][1][0] * entries[1][2][0] == Fraction(1, 4)
    print("weight-one <phi_a|Z_0|phi_b>: diagonal 1/2, off-diagonal -1/2")
    print("<phi_0|P_0 Z_0 P_1 Z_0 P_2|phi_2> = 1/4")

    edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    schedules = list(permutations(edges))
    all_full = [s for s in schedules if all(len(backward_support(s, j)) == 4
                                              for j in range(4))]
    assert len(schedules) == 24 and len(all_full) == 8
    print(f"four-cycle schedules satisfying all-output full support: {len(all_full)}/24")
    print("example surviving schedule:", all_full[0])


if __name__ == "__main__":
    main()
