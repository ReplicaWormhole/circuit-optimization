"""Exact cyclotomic check for four-qubit circuits with rational-pi angles.

This verifies U V4 = D U over a cyclotomic number field. It is independent of
the floating-point checker and accepts any output eigenvalue order.
"""

import argparse
import json
from fractions import Fraction
from math import lcm
from pathlib import Path

import sympy as sp

from check_circuit import ANGLE


def rational_pi(value):
    if not isinstance(value, str):
        raise ValueError("exact checking requires angle strings such as '3*pi/8'")
    match = ANGLE.fullmatch(value.replace(" ", ""))
    if match is None:
        raise ValueError(f"angle is not a rational multiple of pi: {value!r}")
    sign, numerator, denominator = match.groups()
    denominator = int(denominator or 1)
    if denominator == 0:
        raise ValueError("angle denominator must be nonzero")
    return (-1 if sign == "-" else 1) * Fraction(numerator or 1) / denominator


def exact_eigenvalue_labels(candidate):
    if candidate.get("n") != 4:
        raise ValueError("exact checker currently supports only n=4")
    gates = candidate["gates"]
    angle_values = [rational_pi(gate[key]) for gate in gates
                    for key in ("theta", "phi", "lam") if key in gate]
    conductor = lcm(8, *(4 * value.denominator for value in angle_values))
    if sp.totient(conductor) > 128:
        raise ValueError("cyclotomic degree exceeds 128; use a simpler exact ansatz")
    field = sp.QQ.cyclotomic_field(conductor)
    zeta = field.unit
    zero, one = field.zero, field.one
    half = field.from_sympy(sp.Rational(1, 2))
    imaginary = zeta ** (conductor // 4)
    hadamard_value = (zeta ** (conductor // 8) + zeta ** (-conductor // 8)) * half

    def phase(angle, half_angle=False):
        value = rational_pi(angle)
        denominator = (4 if half_angle else 2) * value.denominator
        return zeta ** (value.numerator * (conductor // denominator))

    def local_matrix(gate):
        name = gate["gate"].lower()
        if name == "h":
            h = hadamard_value
            return ((h, h), (h, -h))
        if name == "x":
            return ((zero, one), (one, zero))
        if name == "s":
            return ((one, zero), (zero, imaginary))
        if name == "sdg":
            return ((one, zero), (zero, -imaginary))
        if name in ("rx", "ry", "rz", "u3"):
            positive = phase(gate["theta"], half_angle=True)
            negative = positive ** -1
            cosine = (positive + negative) * half
            sine = (positive - negative) / (2 * imaginary)
            if name == "rx":
                return ((cosine, -imaginary * sine),
                        (-imaginary * sine, cosine))
            if name == "ry":
                return ((cosine, -sine), (sine, cosine))
            if name == "rz":
                return ((negative, zero), (zero, positive))
            phi, lam = phase(gate["phi"]), phase(gate["lam"])
            return ((cosine, -lam * sine),
                    (phi * sine, phi * lam * cosine))
        raise ValueError(f"unsupported gate: {name}")

    unitary = [[one if row == column else zero for column in range(16)]
               for row in range(16)]
    cnot_count = 0
    for position, gate in enumerate(gates):
        name = gate["gate"].lower()
        if name == "cx":
            control, target = gate["control"], gate["target"]
            if (not isinstance(control, int) or not isinstance(target, int)
                    or control == target or not 0 <= control < 4
                    or not 0 <= target < 4):
                raise ValueError(f"gate {position}: invalid CNOT wires")
            control_mask, target_mask = 1 << (3 - control), 1 << (3 - target)
            next_rows = [None] * 16
            for row in range(16):
                output = row ^ target_mask if row & control_mask else row
                next_rows[output] = unitary[row]
            unitary = next_rows
            cnot_count += 1
        else:
            qubit = gate["qubit"]
            if not isinstance(qubit, int) or not 0 <= qubit < 4:
                raise ValueError(f"gate {position}: invalid qubit")
            matrix = local_matrix(gate)
            mask = 1 << (3 - qubit)
            for row in range(16):
                if row & mask:
                    continue
                partner = row ^ mask
                first, second = unitary[row], unitary[partner]
                unitary[row] = [matrix[0][0] * x + matrix[0][1] * y
                                for x, y in zip(first, second)]
                unitary[partner] = [matrix[1][0] * x + matrix[1][1] * y
                                    for x, y in zip(first, second)]

    roots = (one, imaginary, -one, -imaginary)
    labels = []
    for row in unitary:
        matches = [index for index, root in enumerate(roots)
                   if all(row[((basis & 1) << 3) | (basis >> 1)] == root * row[basis]
                          for basis in range(16))]
        if len(matches) != 1:
            return {"exact_diagonalizer": False, "cnot_count": cnot_count,
                    "cyclotomic_conductor": conductor, "output_labels": labels}
        labels.append(matches[0])
    return {"exact_diagonalizer": True, "cnot_count": cnot_count,
            "cyclotomic_conductor": conductor,
            "output_labels": labels,
            "eigenvalue_multiplicities": [labels.count(index) for index in range(4)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    try:
        result = exact_eigenvalue_labels(json.loads(args.candidate.read_text()))
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["exact_diagonalizer"] else 1)


if __name__ == "__main__":
    main()
