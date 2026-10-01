"""Search two reversible quadratic shears for shared V4 orbit coordinates.

One shear cannot make all three length-4 orbits parallel affine planes
(run106). Enumerate distinct shear permutations, then all ordered pairs.
A parallel-plane result permits a common two-bit Fourier position register,
but gives no CNOT count until the nonlinear preprocessor is synthesized.
"""

import json
from itertools import combinations
from pathlib import Path

from algebraic14_orbit_shear import orbits, parity, plane_direction, shear


ROOT = Path(__file__).parent


def form_cost(form):
    v, l1, l2, b1, b2 = form
    return (v.bit_count() + l1.bit_count() + l2.bit_count(),
            b1 + b2, v, l1, l2)


def unique_shears():
    representatives = {}
    for l1, l2 in combinations(range(1, 16), 2):
        for v in range(1, 16):
            if parity(l1 & v) or parity(l2 & v):
                continue
            for b1 in range(2):
                for b2 in range(2):
                    form = (v, l1, l2, b1, b2)
                    permutation = tuple(
                        shear(x, v, l1, l2, b1, b2) for x in range(16))
                    old = representatives.get(permutation)
                    if old is None or form_cost(form) < form_cost(old):
                        representatives[permutation] = form
    return representatives


def main():
    cycles = [cycle for cycle in orbits() if len(cycle) == 4]
    representatives = unique_shears()
    items = list(representatives.items())
    solutions = []
    trials = 0
    all_plane = 0
    parallel_count = 0
    direct_count = 0
    best = None
    for first_perm, first_form in items:
        first_cycles = [[first_perm[x] for x in cycle] for cycle in cycles]
        for second_perm, second_form in items:
            trials += 1
            transformed = [[second_perm[x] for x in cycle]
                           for cycle in first_cycles]
            directions = [plane_direction(cycle)
                          for cycle in transformed]
            if any(direction is None for direction in directions):
                continue
            all_plane += 1
            if not (directions[0] == directions[1] == directions[2]):
                continue
            parallel_count += 1
            cost = form_cost(first_form)[0] + form_cost(second_form)[0]
            direct_tofo = all(
                v.bit_count() == l1.bit_count() == l2.bit_count() == 1
                for v, l1, l2, _, _ in (first_form, second_form))
            direct_count += direct_tofo
            record = {"first": first_form, "second": second_form,
                      "score": cost,
                      "both_direct_toffoli": direct_tofo,
                      "direction": directions[0],
                      "transformed_orbits": transformed}
            if best is None or (cost, record["first"], record["second"]) < (
                    best["score"], best["first"], best["second"]):
                best = record
            if len(solutions) < 16:
                solutions.append(record)
    result = {"unique_shears": len(items),
              "ordered_pairs_tested": trials,
              "all_plane_pairs": all_plane,
              "parallel_plane_pairs": parallel_count,
              "direct_toffoli_pairs": direct_count,
              "best": best,
              "first_solutions": solutions}
    (ROOT / "algebraic14_orbit_two_shears_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
