"""Test whether a diagonal quadratic gauge conjugates fermion translation to V4.

A fermionic cyclic mode shift carries sign
s(x)=x3*(x0+x1+x2) mod 2 when the last occupied mode crosses the rest.
For G|x>=(-1)^Q(x)|x>, equality G Tf G = V4 requires
Q(shift x)+Q(x)=s(x). Enumerate every degree<=2 Boolean polynomial Q.
The orbit holonomy also gives an all-degree obstruction if nonzero.
"""

import json
from itertools import combinations
from pathlib import Path

from algebraic14_orbit_shear import orbits, parity, rotate


ROOT = Path(__file__).parent
MONOMIALS = [(q,) for q in range(4)] + list(combinations(range(4), 2))


def bit(x, q):
    return (x >> (3 - q)) & 1


def monomial(x, wires):
    return all(bit(x, q) for q in wires)


def polynomial(x, coeffs):
    return sum((coeffs >> index) & 1 for index, wires in
               enumerate(MONOMIALS) if monomial(x, wires)) & 1


def fermion_sign(x):
    return bit(x, 3) & parity(x >> 1)


def main():
    cycles = orbits()
    holonomies = [{"orbit": cycle,
                   "signs": [fermion_sign(x) for x in cycle],
                   "holonomy": sum(fermion_sign(x) for x in cycle) & 1}
                  for cycle in cycles]
    solutions = []
    for coeffs in range(1 << len(MONOMIALS)):
        if all(polynomial(rotate(x), coeffs) ^ polynomial(x, coeffs)
               == fermion_sign(x) for x in range(16)):
            solutions.append(coeffs)
    result = {"monomials": MONOMIALS,
              "quadratic_polynomials_tested": 1 << len(MONOMIALS),
              "degree_at_most_two_solutions": len(solutions),
              "orbit_holonomies": holonomies,
              "all_degree_diagonal_gauge_possible": all(
                  row["holonomy"] == 0 for row in holonomies),
              "interpretation":
                  "Odd holonomy on the alternating 2-cycle and 1111 fixed point obstructs any diagonal gauge; fermionic boundary condition must depend on particle parity."}
    (ROOT / "algebraic14_fermion_gauge_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
