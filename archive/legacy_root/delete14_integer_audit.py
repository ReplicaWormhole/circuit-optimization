"""Independent exact integer-polynomial audit for rational-angle V4 circuits.

Arithmetic is in Z[z]/(z^32-z^16+1), where z=exp(i*pi/48) is a primitive
96th root of unity. Every matrix element has one common power-of-two
denominator. Only angles whose half-angle is a multiple of pi/48 are accepted.
The result tests U V4 = D U row by row without SymPy or floating point.
"""

import argparse
from fractions import Fraction
import json
import re


DEGREE = 32
CONDUCTOR = 96
ZERO = (0,) * DEGREE
ONE = (1,) + (0,) * (DEGREE - 1)
ANGLE = re.compile(r"([+-]?)(?:(\d+)\*?)?pi(?:/(\d+))?")


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def multiply_z(a):
    """Multiply by z and reduce z^32 to z^16-1."""
    out = [0] * DEGREE
    for k in range(DEGREE - 1):
        out[k + 1] = a[k]
    out[0] -= a[DEGREE - 1]
    out[16] += a[DEGREE - 1]
    return tuple(out)


POWERS = [ONE]
for _ in range(1, CONDUCTOR):
    POWERS.append(multiply_z(POWERS[-1]))
assert multiply_z(POWERS[-1]) == ONE


def shift(a, power):
    out = [0] * DEGREE
    for j, coefficient in enumerate(a):
        if coefficient:
            for k, monomial_coefficient in enumerate(POWERS[(j + power) % CONDUCTOR]):
                if monomial_coefficient:
                    out[k] += coefficient * monomial_coefficient
    return tuple(out)


def parse_angle(value):
    if not isinstance(value, str):
        raise ValueError(f"expected exact rational-pi string, got {value!r}")
    match = ANGLE.fullmatch(value.replace(" ", ""))
    if match is None:
        raise ValueError(f"unsupported exact angle {value!r}")
    sign, number, denominator = match.groups()
    d = int(denominator or 1)
    if d == 0:
        raise ValueError("zero denominator")
    return Fraction((-1 if sign == "-" else 1) * int(number or 1), d)


def exponent(value, half=False):
    multiple = parse_angle(value) * (24 if half else 48)
    if multiple.denominator != 1:
        raise ValueError(f"{value!r} is outside the zeta96 ring")
    return multiple.numerator


def matrix(gate):
    """Return integer-polynomial 2x2 numerator and denominator exponent."""
    name = gate["gate"].lower()
    if name == "x":
        return ((ZERO, ONE), (ONE, ZERO)), 0
    if name == "s":
        return ((ONE, ZERO), (ZERO, shift(ONE, 24))), 0
    if name == "sdg":
        return ((ONE, ZERO), (ZERO, shift(ONE, 72))), 0
    if name == "h":
        h = add(shift(ONE, 12), shift(ONE, -12))
        return ((h, h), (h, neg(h))), 1
    if name in ("rx", "ry", "rz", "u3"):
        m = exponent(gate["theta"], half=True)
        positive, negative = shift(ONE, m), shift(ONE, -m)
        if name == "rz":
            return ((negative, ZERO), (ZERO, positive)), 0
        cosine = add(positive, negative)
        diff = sub(positive, negative)
        if name == "rx":
            return ((cosine, neg(diff)), (neg(diff), cosine)), 1
        sine = shift(diff, 72)  # divide by i = multiply by -i
        if name == "ry":
            return ((cosine, neg(sine)), (sine, cosine)), 1
        phi, lam = exponent(gate["phi"]), exponent(gate["lam"])
        return ((cosine, neg(shift(sine, lam))),
                (shift(sine, phi), shift(cosine, phi + lam))), 1
    raise ValueError(f"unsupported gate {name!r}")


def check(candidate):
    if candidate.get("n") != 4:
        raise ValueError("requires n=4")
    rows = [[ONE if i == j else ZERO for j in range(16)] for i in range(16)]
    dyadic_power = 0
    cx_count = 0
    for gate in candidate["gates"]:
        if gate["gate"].lower() == "cx":
            control, target = gate["control"], gate["target"]
            if not all(isinstance(q, int) and 0 <= q < 4 for q in (control, target)) or control == target:
                raise ValueError("invalid CNOT wires")
            cx_count += 1
            control_mask, target_mask = 8 >> control, 8 >> target
            for row in range(16):
                if row & control_mask and not row & target_mask:
                    rows[row], rows[row ^ target_mask] = rows[row ^ target_mask], rows[row]
            continue
        q = gate["qubit"]
        if not isinstance(q, int) or not 0 <= q < 4:
            raise ValueError("invalid local wire")
        local, scale = matrix(gate)
        (m00, m01), (m10, m11) = local
        mask = 8 >> q
        for row in range(16):
            if row & mask:
                continue
            partner = row ^ mask
            upper, lower = rows[row], rows[partner]
            next_upper = []
            next_lower = []
            for a, b in zip(upper, lower):
                # Gate entries have at most two z-monomials. Multiplication by
                # these fixed entries uses sparse integer shifts only.
                def product(poly, value):
                    result = ZERO
                    for k, coefficient in enumerate(poly):
                        if coefficient:
                            term = shift(value, k)
                            result = add(result, tuple(coefficient * x for x in term))
                    return result
                next_upper.append(add(product(m00, a), product(m01, b)))
                next_lower.append(add(product(m10, a), product(m11, b)))
            rows[row], rows[partner] = next_upper, next_lower
        dyadic_power += scale
    labels = []
    for row in rows:
        matches = []
        for name, phase_power in (("1", 0), ("i", 24), ("-1", 48), ("-i", 72)):
            if all(row[((col & 1) << 3) | (col >> 1)] == shift(row[col], phase_power)
                   for col in range(16)):
                matches.append(name)
        if len(matches) != 1:
            return {"exact_uv_equals_du": False, "cnot_count": cx_count,
                    "failed_row": len(labels), "matches": matches,
                    "dyadic_denominator_power": dyadic_power}
        labels.append(matches[0])
    return {"exact_uv_equals_du": True, "cnot_count": cx_count,
            "cyclotomic_conductor": CONDUCTOR,
            "dyadic_denominator_power": dyadic_power,
            "eigenvalue_labels": labels}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate")
    args = parser.parse_args()
    print(json.dumps(check(json.load(open(args.candidate))), indent=2))
