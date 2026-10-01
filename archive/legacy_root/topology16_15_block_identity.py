"""Exact Q(zeta_32) proofs of both two-CNOT butterfly replacements."""

import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent
F = sp.QQ.cyclotomic_field(32)
Z = F.unit
ZERO, ONE = F.zero, F.one
HALF = F.from_sympy(sp.Rational(1, 2))
I = Z ** 8


def local(gate):
    angle = rational_pi(gate["theta"])
    assert 8 % angle.denominator == 0
    p = Z ** (angle.numerator * (8 // angle.denominator))
    c = (p + p ** -1) * HALF
    s = (p - p ** -1) * HALF / I
    kind = gate["gate"]
    if kind == "rx":
        return ((c, -I * s), (-I * s, c))
    if kind == "ry":
        return ((c, -s), (s, c))
    if kind == "rz":
        return ((p ** -1, ZERO), (ZERO, p))
    raise ValueError(gate)


def gate_matrix(gate, upper, lower):
    result = [[ZERO for _ in range(4)] for _ in range(4)]
    if gate["gate"] == "cx":
        control = 0 if gate["control"] == upper else 1
        target = 0 if gate["target"] == upper else 1
        for x in range(4):
            y = x ^ (2 if target == 0 else 1) if x & (2 if control == 0 else 1) else x
            result[y][x] = ONE
        return result
    unit = local(gate)
    axis = 0 if gate["qubit"] == upper else 1
    mask = 2 if axis == 0 else 1
    for x in range(4):
        for output_bit in range(2):
            y = (x & ~mask) | (output_bit * mask)
            input_bit = 1 if x & mask else 0
            result[y][x] = unit[output_bit][input_bit]
    return result


def multiply(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(4)), ZERO)
             for j in range(4)] for i in range(4)]


def circuit_matrix(gates, upper, lower):
    result = [[ONE if i == j else ZERO for j in range(4)] for i in range(4)]
    for gate in gates:
        result = multiply(gate_matrix(gate, upper, lower), result)
    return result


def equal_up_to_known_phase(source, replacement, phase, upper, lower):
    left = circuit_matrix(source, upper, lower)
    right = circuit_matrix(replacement, upper, lower)
    return all(right[i][j] == phase * left[i][j]
               for i in range(4) for j in range(4))


def main():
    gates = json.loads((ROOT / "baseline_18.json").read_text())["gates"]
    source_forward = [{"gate": "cx", "control": 0, "target": 2}, *gates[18:28]]
    source_reverse = [{"gate": "cx", "control": 3, "target": 1}, *gates[28:38]]
    forward = equal_up_to_known_phase(source_forward, exact_block(), Z ** 12, 0, 2)
    reverse = equal_up_to_known_phase(source_reverse, reverse_block(), Z ** -4, 1, 3)
    assert forward and reverse
    print(json.dumps({"field": "Q(zeta_32)",
                      "forward_identity": forward,
                      "forward_phase": "exp(3*pi*i/4)",
                      "reverse_identity": reverse,
                      "reverse_phase": "exp(-pi*i/4)"}))


if __name__ == "__main__":
    main()
