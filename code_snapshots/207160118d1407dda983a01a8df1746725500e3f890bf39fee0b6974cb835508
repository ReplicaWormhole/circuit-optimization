"""Exact cut-rank screen for one output-CNOT gauge on the terminal suffix.

S = F_01 tensor F_23, F = exp(i*pi*(XX+YY)/8), has four-CNOT synthesis.
Every output CNOT is a computational-basis permutation, so P*S is still a
valid suffix gauge. Check all 12 directed CNOTs to see whether any P*S can
use at most three CNOTs. A rank above eight across any bipartition rules
that out because each crossing CNOT multiplies Schmidt rank by at most two.
This does not address gauges spanning earlier circuit boundaries.
"""

import json
import math

import sympy as s


I = s.I
R2 = s.sqrt(2)
F = s.Matrix([[1, 0, 0, 0],
              [0, 1/R2, I/R2, 0],
              [0, I/R2, 1/R2, 0],
              [0, 0, 0, 1]])
S = s.kronecker_product(F, F)
CUTS = ((0, 1), (0, 2), (0, 3))


def bit(x, wire):
    return (x >> (3 - wire)) & 1


def cnot(control, target):
    matrix = s.zeros(16)
    for x in range(16):
        y = x ^ ((1 << (3 - target)) if bit(x, control) else 0)
        matrix[y, x] = 1
    return matrix


def realign(unitary, a):
    b = tuple(w for w in range(4) if w not in a)
    result = s.zeros(16)
    for output in range(16):
        for input_ in range(16):
            oa = 2 * bit(output, a[0]) + bit(output, a[1])
            ia = 2 * bit(input_, a[0]) + bit(input_, a[1])
            ob = 2 * bit(output, b[0]) + bit(output, b[1])
            ib = 2 * bit(input_, b[0]) + bit(input_, b[1])
            result[4 * oa + ia, 4 * ob + ib] = unitary[output, input_]
    return result


def main():
    rows = []
    for control in range(4):
        for target in range(4):
            if control == target:
                continue
            unitary = cnot(control, target) * S
            ranks = {"".join(map(str, a)): int(realign(unitary, a).rank())
                     for a in CUTS}
            required = max(math.ceil(math.log2(rank)) for rank in ranks.values())
            rows.append({"output_cnot": [control, target],
                         "cross_pair": (control < 2) != (target < 2),
                         "cut_ranks": ranks,
                         "three_cnot_suffix_excluded": required > 3})
    print(json.dumps({"suffix": "F01(pi/8,pi/8) F23(pi/8,pi/8)",
                      "total_directed_output_cnots": len(rows),
                      "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
