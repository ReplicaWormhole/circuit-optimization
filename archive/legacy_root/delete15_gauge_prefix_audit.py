"""Independent exact phase audit of the six-CNOT gauge-1 prefix.

Tracks computational-basis states and exponents of z=exp(i*pi/16).
Checks the prefix's endpoint map and diagonal parity polynomial, and checks
that the added triple-orbit gauge commutes with the right cyclic shift.
"""

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANGLE = re.compile(r"([+-]?)(?:(\d+)\*?)?pi(?:/(\d+))?")
TRIPLES = (0b1011, 0b1101, 0b1110, 0b0111)


def eighths(value):
    match = ANGLE.fullmatch(value)
    if not match:
        raise ValueError(f"unsupported angle {value!r}")
    sign, factor, divisor = match.groups()
    number = (-1 if sign == "-" else 1) * int(factor or 1) * 8
    divisor = int(divisor or 1)
    if number % divisor:
        raise ValueError(f"angle not a multiple of pi/8: {value!r}")
    return number // divisor


def parity_sign(x, mask):
    return -1 if (x & mask).bit_count() % 2 else 1


def right_shift(x):
    return ((x & 1) << 3) | (x >> 1)


def linear_rows(gates):
    rows = [8, 4, 2, 1]
    for gate in gates:
        if gate["gate"] == "cx":
            rows[gate["target"]] ^= rows[gate["control"]]
    return rows


def actual_phase_and_output(gates, x):
    state = x
    power = 0
    for gate in gates:
        if gate["gate"] == "cx":
            c, t = 8 >> gate["control"], 8 >> gate["target"]
            if state & c:
                state ^= t
        elif gate["gate"] == "rz":
            bit = 8 >> gate["qubit"]
            power += eighths(gate["theta"]) * (1 if state & bit else -1)
        else:
            raise ValueError(f"unexpected gate in prefix: {gate['gate']}")
    return state, power % 32


def target_phase(x):
    terms = [(0b0100, 5), (0b0010, 6), (0b0001, 3),
             (0b1011, -1), (0b1110, 1), (0b0111, -2), (0b0110, -4)]
    return sum(-coefficient * parity_sign(x, mask)
               for mask, coefficient in terms) % 32


def gauge_phase(x):
    # Relative to the 17-CNOT target (triple angles 1,2,3,0 eighths),
    # gauge 1 subtracts 2 eighths on each of the four triple masks.
    return sum(2 * parity_sign(x, mask) for mask in TRIPLES) % 32


def main():
    candidate = json.loads((ROOT / "algebraic16_to15_gauge_1_best.json").read_text())
    prefix = candidate["gates"][:13]
    if sum(g["gate"] == "cx" for g in prefix) != 6:
        raise AssertionError("expected a six-CNOT prefix")
    rows = linear_rows(prefix)
    offsets = set()
    for x in range(16):
        actual_output, actual_power = actual_phase_and_output(prefix, x)
        expected_output = sum((((x & mask).bit_count() % 2) << (3 - q))
                              for q, mask in enumerate(rows))
        if actual_output != expected_output:
            raise AssertionError(f"basis {x}: wrong endpoint map")
        offsets.add((actual_power - target_phase(x)) % 32)
        if gauge_phase(x) != gauge_phase(right_shift(x)):
            raise AssertionError(f"gauge does not commute with shift on basis {x}")
    if len(offsets) != 1:
        raise AssertionError(f"non-global prefix phase mismatch: {offsets}")
    print(json.dumps({"prefix_cnot_count": 6, "endpoint_rows": rows,
                      "exact_phase_match_up_to_global": True,
                      "global_phase_exponent_mod32": next(iter(offsets)),
                      "triple_orbit_gauge_commutes_with_right_shift": True}, indent=2))


if __name__ == "__main__":
    main()
