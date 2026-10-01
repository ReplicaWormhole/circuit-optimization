"""Induced V4 action after one quadratic shear and a two-CNOT linear map.

Use run106's representative S=(15,3,5,0,0), making two length-four orbits
parallel affine planes. Then apply CX0->2,CX1->2. The conserved parity bit
is z=y2 xor y3. Derive exact 8-state actions in the two z sectors.
"""

import argparse
import json
from pathlib import Path

from algebraic14_orbit_shear import orbits, rotate, shear
from topology14_orbit_action import linear_map


ROOT = Path(__file__).resolve().parent
SHEAR = (15, 3, 5, 0, 0)


def parity_sector(y):
    return ((y >> 1) ^ y) & 1


def sector_index(y):
    return y >> 1


def cycles(permutation):
    remaining = set(range(len(permutation)))
    result = []
    while remaining:
        initial = min(remaining)
        orbit = [initial]
        current = permutation[initial]
        while current != initial:
            orbit.append(current)
            current = permutation[current]
        remaining.difference_update(orbit)
        result.append(orbit)
    return result


def anf(truth):
    coeff = truth[:]
    n = len(truth).bit_length() - 1
    for bit in range(n):
        for mask in range(1 << n):
            if mask & (1 << bit):
                coeff[mask] ^= coeff[mask ^ (1 << bit)]
    return [mask for mask, value in enumerate(coeff) if value]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    coordinate = [linear_map(shear(x, *SHEAR)) for x in range(16)]
    assert len(set(coordinate)) == 16
    inverse = [coordinate.index(y) for y in range(16)]
    induced = [coordinate[rotate(inverse[y])] for y in range(16)]
    assert all(parity_sector(y) == parity_sector(induced[y]) for y in range(16))
    sectors = {}
    for z in (0, 1):
        domain = [y for y in range(16) if parity_sector(y) == z]
        assert [sector_index(y) for y in domain] == list(range(8))
        permutation = [sector_index(induced[y]) for y in domain]
        sectors[str(z)] = {"permutation": permutation, "cycles": cycles(permutation),
                           "output_bit_anf_masks": [
                               anf([(value >> (2 - bit)) & 1 for value in permutation])
                               for bit in range(3)]}
    result = {"shear": SHEAR,
              "linear_postmap": "CX0->2,CX1->2",
              "coordinate_map": coordinate,
              "induced_permutation": induced,
              "global_parity_equals_y2_xor_y3": True,
              "transformed_length4_orbits": [[coordinate[x] for x in cycle]
                                             for cycle in orbits() if len(cycle) == 4],
              "sectors": sectors}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
