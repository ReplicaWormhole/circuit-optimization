"""Independent exact arithmetic audit for four-qubit pi/8-angle circuits.

Uses the cyclotomic ring Z[z]/(z^16+1), z=exp(i*pi/16), with a global dyadic
denominator. Checks rowwise U V = D U for D in {1,i,-1,-i}.
"""

import argparse
import json
import re


ZERO = (0,) * 16
ONE = (1,) + (0,) * 15
RX_PATTERN = re.compile(r"([+-]?)(?:(\d+)\*?)?pi(?:/(\d+))?")


def shift(v, k):
    result = [0] * 16
    for j, coefficient in enumerate(v):
        power = (j + k) % 32
        if power >= 16:
            result[power - 16] -= coefficient
        else:
            result[power] += coefficient
    return tuple(result)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def exponent(angle, denominator):
    match = RX_PATTERN.fullmatch(str(angle))
    if not match:
        raise ValueError(f"unsupported exact angle {angle!r}")
    sign, multiplier, div = match.groups()
    numerator = (-1 if sign == "-" else 1) * int(multiplier or 1) * denominator
    div = int(div or 1)
    if numerator % div:
        raise ValueError(f"angle {angle!r} is not a multiple of pi/{denominator}")
    return numerator // div


def check(candidate):
    if candidate["n"] != 4:
        raise ValueError("requires four qubits")
    # state[row][column] is the exact matrix entry numerator.
    state = [[ONE if row == column else ZERO for column in range(16)]
             for row in range(16)]
    denominator_power = 0
    cx_count = 0
    for gate in candidate["gates"]:
        name = gate["gate"]
        if name == "cx":
            cx_count += 1
            cmask = 8 >> gate["control"]
            tmask = 8 >> gate["target"]
            for row in range(16):
                if row & cmask and not row & tmask:
                    state[row], state[row ^ tmask] = state[row ^ tmask], state[row]
        elif name == "rz":
            qmask = 8 >> gate["qubit"]
            m = exponent(gate["theta"], 8)
            for row in range(16):
                state[row] = [shift(value, m if row & qmask else -m)
                              for value in state[row]]
        elif name in ("rx", "ry", "h"):
            qmask = 8 >> gate["qubit"]
            k = exponent(gate["theta"], 4) if name in ("rx", "ry") else None
            for row in range(16):
                if row & qmask:
                    continue
                upper, lower = state[row], state[row ^ qmask]
                next_upper, next_lower = [], []
                for a, b in zip(upper, lower):
                    if name in ("rx", "ry"):
                        ca = add(shift(a, 2*k), shift(a, -2*k))
                        cb = add(shift(b, 2*k), shift(b, -2*k))
                        if name == "rx":
                            oa = sub(shift(a, -2*k), shift(a, 2*k))
                            ob = sub(shift(b, -2*k), shift(b, 2*k))
                            next_upper.append(add(ca, ob))
                            next_lower.append(add(oa, cb))
                        else:
                            sa = shift(sub(shift(a, 2*k), shift(a, -2*k)), 24)
                            sb = shift(sub(shift(b, 2*k), shift(b, -2*k)), 24)
                            next_upper.append(sub(ca, sb))
                            next_lower.append(add(sa, cb))
                    else:
                        ha = add(shift(a, 4), shift(a, -4))
                        hb = add(shift(b, 4), shift(b, -4))
                        next_upper.append(add(ha, hb))
                        next_lower.append(sub(ha, hb))
                state[row], state[row ^ qmask] = next_upper, next_lower
            denominator_power += 1
        else:
            raise ValueError(f"unsupported gate for exact audit: {name}")
    labels = []
    for row in range(16):
        possible = []
        for label, power in (("1", 0), ("i", 8), ("-1", 16), ("-i", 24)):
            if all(state[row][((col & 1) << 3) | (col >> 1)] == shift(state[row][col], power)
                   for col in range(16)):
                possible.append(label)
        if len(possible) != 1:
            raise AssertionError(f"row {row} has {len(possible)} eigenvalue labels: {possible}")
        labels.append(possible[0])
    return {"cnot_count": cx_count, "exact_uv_equals_du": True,
            "dyadic_denominator_power": denominator_power, "eigenvalue_labels": labels}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate")
    args = parser.parse_args()
    print(json.dumps(check(json.load(open(args.candidate))), indent=2))
