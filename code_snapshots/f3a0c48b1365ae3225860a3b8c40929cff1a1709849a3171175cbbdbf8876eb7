"""Independently check the exact 15-CNOT gate list over Q(zeta_32)."""

import argparse
import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi


ROOT = Path(__file__).parent
MATRIX_STRINGS = [["1/2", "(sqrt(2)+i)/2"],
                  ["(-sqrt(2)+i)/2", "1/2"]]


def check(spec):
    assert spec["n"] == 4
    field = sp.QQ.cyclotomic_field(32)
    z = field.unit
    zero, one = field.zero, field.one
    half = field.from_sympy(sp.Rational(1, 2))
    imaginary = z ** 8
    sqrt2 = z ** 4 + z ** -4

    def phase(value, half_angle=False):
        coefficient = rational_pi(value)
        divisor = (4 if half_angle else 2) * coefficient.denominator
        if 32 % divisor:
            raise ValueError(f"angle outside Q(zeta_32): {value}")
        return z ** (coefficient.numerator * (32 // divisor))

    def local(gate):
        name = gate["gate"]
        if name == "u3_algebraic":
            assert gate["qubit"] in (0, 1)
            assert gate["theta"] == "2*pi/3"
            assert gate["phi"] == "pi-atan(1/sqrt(2))"
            assert gate["lam"] == "-pi+atan(1/sqrt(2))"
            assert gate["matrix"] == MATRIX_STRINGS
            return ((half, (sqrt2 + imaginary) * half),
                    ((-sqrt2 + imaginary) * half, half))
        if name == "h":
            h = sqrt2 * half
            return ((h, h), (h, -h))
        if name == "x":
            return ((zero, one), (one, zero))
        if name == "s":
            return ((one, zero), (zero, imaginary))
        if name == "sdg":
            return ((one, zero), (zero, -imaginary))
        if name in ("rx", "ry", "rz", "u3"):
            p = phase(gate["theta"], half_angle=True)
            q = p ** -1
            cosine = (p + q) * half
            sine = (p - q) / (2 * imaginary)
            if name == "rx":
                return ((cosine, -imaginary * sine),
                        (-imaginary * sine, cosine))
            if name == "ry":
                return ((cosine, -sine), (sine, cosine))
            if name == "rz":
                return ((q, zero), (zero, p))
            phi, lam = phase(gate["phi"]), phase(gate["lam"])
            return ((cosine, -lam * sine),
                    (phi * sine, phi * lam * cosine))
        raise ValueError(name)

    rows = [[one if i == j else zero for j in range(16)] for i in range(16)]
    count = 0
    for position, gate in enumerate(spec["gates"]):
        if gate["gate"] == "cx":
            control, target = gate["control"], gate["target"]
            assert 0 <= control < 4 and 0 <= target < 4 and control != target
            cmask, tmask = 1 << (3 - control), 1 << (3 - target)
            next_rows = [None] * 16
            for row in range(16):
                output = row ^ tmask if row & cmask else row
                next_rows[output] = rows[row]
            rows = next_rows
            count += 1
        else:
            qubit = gate["qubit"]
            assert 0 <= qubit < 4
            matrix = local(gate)
            mask = 1 << (3 - qubit)
            for row in range(16):
                if row & mask:
                    continue
                partner = row ^ mask
                first, second = rows[row], rows[partner]
                rows[row] = [matrix[0][0] * a + matrix[0][1] * b
                             for a, b in zip(first, second)]
                rows[partner] = [matrix[1][0] * a + matrix[1][1] * b
                                 for a, b in zip(first, second)]
    roots = (one, imaginary, -one, -imaginary)
    labels = []
    for row in rows:
        matches = [index for index, root in enumerate(roots)
                   if all(row[((basis & 1) << 3) | (basis >> 1)] == root * row[basis]
                          for basis in range(16))]
        if len(matches) != 1:
            return {"exact_diagonalizer": False, "cnot_count": count,
                    "output_labels_so_far": labels}
        labels.append(matches[0])
    return {"exact_diagonalizer": True, "cnot_count": count,
            "cyclotomic_conductor": 32, "output_labels": labels,
            "eigenvalue_multiplicities": [labels.count(k) for k in range(4)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", nargs="?", type=Path,
                        default=ROOT / "topology16_15_exact_candidate.json")
    args = parser.parse_args()
    spec = json.loads(args.candidate.read_text())
    result = check(spec)
    print(json.dumps(result, indent=2))
    if not result["exact_diagonalizer"] or result["cnot_count"] != 15:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
