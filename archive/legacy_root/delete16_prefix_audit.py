"""Independent finite-integer audit of the seven-CNOT parity prefix.

All Rz angles are multiples of pi/8. Their phases are powers of
z=exp(i*pi/16), so integer exponents modulo 32 suffice. The expected
diagonal target is the four parity rotations stated in the 17-CNOT
derivation, followed by the endpoint CNOT 0->2.
"""

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANGLE = re.compile(r"([+-]?)(?:(\d+)\*?)?pi(?:/(\d+))?")


def eighths(value):
    match = ANGLE.fullmatch(value)
    if not match:
        raise ValueError(f"not a rational pi angle: {value!r}")
    sign, factor, divisor = match.groups()
    number = (-1 if sign == "-" else 1) * int(factor or 1) * 8
    divisor = int(divisor or 1)
    if number % divisor:
        raise ValueError(f"not a multiple of pi/8: {value!r}")
    return number // divisor


def run_prefix(gates, x):
    state = x
    power = 0
    for gate in gates:
        if gate["gate"] == "cx":
            control = 8 >> gate["control"]
            target = 8 >> gate["target"]
            if state & control:
                state ^= target
        elif gate["gate"] == "rz":
            target = 8 >> gate["qubit"]
            power += eighths(gate["theta"]) * (1 if state & target else -1)
        else:
            raise ValueError(f"unexpected prefix gate {gate['gate']}")
    return state, power % 32


def target_power(x):
    local = [(0b0100, 5), (0b0010, 6), (0b0001, 3)]
    nonlocal_terms = [(0b1011, 1), (0b1101, 2),
                      (0b1110, 3), (0b0110, -4)]
    return sum(k * (1 if (x & mask).bit_count() % 2 else -1)
               for mask, k in local + nonlocal_terms) % 32


def main():
    exact = json.loads((ROOT / "algebraic16_exact_candidate.json").read_text())
    prefix = exact["gates"][:14]
    if sum(g["gate"] == "cx" for g in prefix) != 7:
        raise AssertionError("prefix CNOT count changed")
    results = []
    offsets = set()
    for x in range(16):
        out, actual_power = run_prefix(prefix, x)
        expected_out = x ^ 0b0010 if x & 0b1000 else x
        if out != expected_out:
            raise AssertionError(f"basis {x}: endpoint maps to {out}, expected {expected_out}")
        desired_power = target_power(x)
        offsets.add((actual_power - desired_power) % 32)
        results.append([x, out, actual_power, desired_power])
    if len(offsets) != 1:
        raise AssertionError(f"prefix phases differ by non-global offsets {offsets}")
    print(json.dumps({"prefix_cnot_count": 7, "endpoint": [0, 2],
                      "exact_mask_and_phase_match_up_to_global": True,
                      "global_phase_exponent_mod_32": next(iter(offsets)),
                      "rows_input_output_actual_desired": results}, indent=2))


if __name__ == "__main__":
    main()
