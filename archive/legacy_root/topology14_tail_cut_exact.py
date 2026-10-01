"""Exact cyclotomic cut-rank certificate for the fixed final Fourier suffix.

Up to input-side one-qubit rotations, the suffix is (F01 tensor F23) CZ12.
Both F blocks are computed from baseline_18 in Q(zeta_16). Exact Gaussian
elimination determines full realignment rank on cuts 02|13 and 03|12.
"""

import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi


ROOT = Path(__file__).resolve().parent
FIELD = sp.QQ.cyclotomic_field(16)
Z = FIELD.unit
ZERO, ONE = FIELD.zero, FIELD.one
HALF = FIELD.from_sympy(sp.Rational(1, 2))
I = Z ** 4


def local(gate):
    value = rational_pi(gate["theta"])
    assert 8 % value.denominator == 0
    p = Z ** (value.numerator * (8 // value.denominator))
    c = (p + p ** -1) * HALF
    s = (p - p ** -1) * HALF / I
    if gate["gate"] == "rx":
        return ((c, -I * s), (-I * s, c))
    if gate["gate"] == "rz":
        return ((p ** -1, ZERO), (ZERO, p))
    raise ValueError(gate)


def two_qubit(gates, first_wire, second_wire):
    u = [[ONE if a == b else ZERO for b in range(4)] for a in range(4)]
    for gate in gates:
        g = [[ZERO] * 4 for _ in range(4)]
        if gate["gate"] == "cx":
            control = 0 if gate["control"] == first_wire else 1
            target = 0 if gate["target"] == first_wire else 1
            for x in range(4):
                y = x ^ (2 if target == 0 else 1) if x & (2 if control == 0 else 1) else x
                g[y][x] = ONE
        else:
            axis = 0 if gate["qubit"] == first_wire else 1
            mask = 2 if axis == 0 else 1
            m = local(gate)
            for x in range(4):
                for output_bit in range(2):
                    y = (x & ~mask) | (output_bit * mask)
                    input_bit = int(bool(x & mask))
                    g[y][x] = m[output_bit][input_bit]
        u = [[sum((g[a][k] * u[k][b] for k in range(4)), ZERO)
              for b in range(4)] for a in range(4)]
    return u


def effective_suffix(f, g):
    result = [[ZERO] * 16 for _ in range(16)]
    for output in range(16):
        o01, o23 = output >> 2, output & 3
        for inp in range(16):
            i01, i23 = inp >> 2, inp & 3
            phase = -ONE if ((inp >> 2) & 1) and ((inp >> 1) & 1) else ONE
            result[output][inp] = f[o01][i01] * g[o23][i23] * phase
    return result


def bits(index, wires):
    return sum(((index >> (3 - wire)) & 1) << (len(wires) - 1 - j)
               for j, wire in enumerate(wires))


def realign(u, left):
    right = tuple(q for q in range(4) if q not in left)
    out = [[ZERO] * 16 for _ in range(16)]
    for output in range(16):
        for inp in range(16):
            row = bits(output, left) * 4 + bits(inp, left)
            col = bits(output, right) * 4 + bits(inp, right)
            out[row][col] = u[output][inp]
    return out


def determinant(matrix):
    matrix = [row[:] for row in matrix]
    answer = ONE
    for pivot in range(16):
        row = next((row for row in range(pivot, 16)
                    if matrix[row][pivot] != ZERO), None)
        if row is None:
            return ZERO
        if row != pivot:
            matrix[pivot], matrix[row] = matrix[row], matrix[pivot]
            answer = -answer
        value = matrix[pivot][pivot]
        answer *= value
        reciprocal = value ** -1
        for row in range(pivot + 1, 16):
            factor = matrix[row][pivot] * reciprocal
            if factor == ZERO:
                continue
            for col in range(pivot + 1, 16):
                matrix[row][col] -= factor * matrix[pivot][col]
            matrix[row][pivot] = ZERO
    return answer


def main():
    gates = json.loads((ROOT / "baseline_18.json").read_text())["gates"]
    f = two_qubit(gates[43:52], 0, 1)
    g = two_qubit(gates[53:62], 2, 3)
    u = effective_suffix(f, g)
    determinants = {}
    for cut in ((0, 2), (0, 3)):
        value = determinant(realign(u, cut))
        assert value != ZERO
        determinants["".join(map(str, cut))] = str(FIELD.to_sympy(value))
    print(json.dumps({"field": "Q(zeta_16)",
                      "effective_operator": "(F01 tensor F23) CZ12; omitted input-side one-qubit gates",
                      "realignment_determinants": determinants,
                      "rank_01_upper_bound": 2,
                      "conclusion": "fixed suffix requires at least five CNOTs"},
                     indent=2))


if __name__ == "__main__":
    main()
