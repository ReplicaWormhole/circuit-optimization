"""Derive V4 action after algebraic run108's best two-shear coordinate map.

The shears make four parallel cosets of H={0,6,10,12}. Apply a two-CNOT
linear map y2=x0 xor x1 xor x2, retaining y0=x0,y1=x1,y3=x3. Thus
position bits are (y0,y1), label bits are (y2,y3). For each label, fit the
induced right-shift action on position bits to an affine F2^2 map.
"""

import argparse
import json
from pathlib import Path

from algebraic14_orbit_shear import rotate, shear


ROOT = Path(__file__).resolve().parent


def apply_shears(x, first, second):
    return shear(shear(x, *first), *second)


def linear_map(x):
    parity012 = ((x >> 3) ^ (x >> 2) ^ (x >> 1)) & 1
    return (x & 0b1101) | (parity012 << 1)


def position(y):
    return ((y >> 3) & 1) * 2 + ((y >> 2) & 1)


def label(y):
    return ((y >> 1) & 1) * 2 + (y & 1)


def affine_fit(mapping):
    shift = mapping[0]
    col0, col1 = mapping[2] ^ shift, mapping[1] ^ shift
    return {"A_columns": [col0, col1], "b": shift,
            "valid": mapping[3] == col0 ^ col1 ^ shift}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = json.loads((ROOT / "algebraic14_orbit_two_shears_result.json").read_text())["best"]
    first, second = record["first"], record["second"]
    coordinate_map = [linear_map(apply_shears(x, first, second)) for x in range(16)]
    assert len(set(coordinate_map)) == 16
    inverse = [coordinate_map.index(y) for y in range(16)]
    conjugated = [coordinate_map[rotate(inverse[y])] for y in range(16)]
    by_label = {}
    for ell in range(4):
        mapping = [None] * 4
        for y, out in enumerate(conjugated):
            if label(y) != ell:
                continue
            assert label(out) == ell
            mapping[position(y)] = position(out)
        by_label[str(ell)] = {"position_permutation": mapping,
                              "affine": affine_fit(mapping)}
    result = {"shears": [first, second],
              "linear_postmap": "CX(0,2), CX(1,2)",
              "coordinate_map": coordinate_map,
              "conjugated_permutation": conjugated,
              "by_label": by_label,
              "all_affine": all(row["affine"]["valid"] for row in by_label.values()),
              "same_linear_part": len({tuple(row["affine"]["A_columns"])
                                        for row in by_label.values()}) == 1}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
