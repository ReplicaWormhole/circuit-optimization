"""Enumerate one quadratic reversible shear as an orbit-coordinate preprocessor.

Any reversible quadratic shear conjugate to a Toffoli by affine maps has form
x -> x xor v (l1(x) xor b1)(l2(x) xor b2), with l1(v)=l2(v)=0.
A single shared 2-bit Fourier position register requires the three size-4
right-shift orbits to map to parallel affine 2-planes. Test that necessary
condition exhaustively for this class. A negative result is only restricted
to one shear and this orbit encoding.
"""

import json
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).parent


def rotate(x):
    return ((x & 1) << 3) | (x >> 1)


def orbits():
    unseen = set(range(16))
    result = []
    while unseen:
        start = min(unseen)
        orbit = [start]
        point = rotate(start)
        while point != start:
            orbit.append(point)
            point = rotate(point)
        unseen.difference_update(orbit)
        result.append(orbit)
    return result


def parity(x):
    return x.bit_count() & 1


def shear(x, v, l1, l2, b1, b2):
    return x ^ (v if (parity(x & l1) ^ b1) and
                (parity(x & l2) ^ b2) else 0)


def plane_direction(orbit):
    if len(orbit) != 4 or orbit[0] ^ orbit[1] ^ orbit[2] ^ orbit[3]:
        return None
    anchor = orbit[0]
    return tuple(sorted(anchor ^ point for point in orbit))


def main():
    cycles = [orbit for orbit in orbits() if len(orbit) == 4]
    assert len(cycles) == 3
    original = [plane_direction(cycle) for cycle in cycles]
    trials = 0
    plane_counts = {str(index): 0 for index in range(4)}
    parallel_examples = []
    all_plane_examples = []
    for l1, l2 in combinations(range(1, 16), 2):
        for v in range(1, 16):
            if parity(l1 & v) or parity(l2 & v):
                continue
            for b1 in range(2):
                for b2 in range(2):
                    trials += 1
                    transformed = [[shear(x, v, l1, l2, b1, b2)
                                    for x in cycle] for cycle in cycles]
                    directions = [plane_direction(cycle)
                                  for cycle in transformed]
                    count = sum(direction is not None
                                for direction in directions)
                    plane_counts[str(count)] += 1
                    if count == 3:
                        example = {"v": v, "l1": l1, "l2": l2,
                                   "b1": b1, "b2": b2,
                                   "directions": directions,
                                   "transformed_orbits": transformed}
                        if len(all_plane_examples) < 4:
                            all_plane_examples.append(example)
                        if directions[0] == directions[1] == directions[2]:
                            parallel_examples.append(example)
    result = {"right_shift_orbits": orbits(),
              "original_plane_directions": original,
              "shears_tested": trials,
              "plane_count_histogram": plane_counts,
              "all_three_plane_examples": all_plane_examples,
              "parallel_shear_count": len(parallel_examples),
              "parallel_examples": parallel_examples[:8]}
    (ROOT / "algebraic14_orbit_shear_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("shears_tested", "plane_count_histogram",
                       "parallel_shear_count", "parallel_examples")},
                     indent=2))


if __name__ == "__main__":
    main()
