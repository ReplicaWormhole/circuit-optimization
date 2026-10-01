"""Symbolic determinant of gauged-boundary realignment across 01|23.

Let t=e^{-i delta/2}. The boundary is (F13 CX3,1)(F02 CX0,2)
RZZ12(delta), up to an irrelevant global phase. The two exact 4x4
butterfly matrices are taken from topology16_15_block_identity.py.
The 16x16 realignment determinant being nonzero proves Schmidt rank 16.
"""

import json

import sympy as sp


S = sp.sqrt(2) / 2
FORWARD = sp.Matrix([
    [-1, 0, 0, 0],
    [0, S, 0, -S],
    [0, -S, 0, -S],
    [0, 0, 1, 0],
])
REVERSE = sp.Matrix([
    [-1, 0, 0, 0],
    [0, 0, -S, S],
    [0, 0, -S, -S],
    [0, 1, 0, 0],
])


def bit(x, q):
    return (x >> (3 - q)) & 1


def realignment(t):
    m = sp.zeros(16)
    for output in range(16):
        o0, o1, o2, o3 = (bit(output, q) for q in range(4))
        for inp in range(16):
            i0, i1, i2, i3 = (bit(inp, q) for q in range(4))
            row = ((o0 * 2 + o1) * 4) + (i0 * 2 + i1)
            col = ((o2 * 2 + o3) * 4) + (i2 * 2 + i3)
            m[row, col] = (FORWARD[o0 * 2 + o2, i0 * 2 + i2]
                           * REVERSE[o1 * 2 + o3, i1 * 2 + i3]
                           * (t if i1 == i2 else 1 / t))
    return m


def main():
    t = sp.symbols("t", nonzero=True)
    m = realignment(t)
    determinant = sp.factor(m.det(method="domain-ge"))
    assert determinant != 0
    print(json.dumps({"nonzero_entries": sum(value != 0 for value in m),
                      "determinant": str(determinant)}))


if __name__ == "__main__":
    main()
