"""Exact balanced-cut determinants for two-body cross-pair parity phases.

Compute both 16x16 realignment determinants of G(t)(F01 tensor F23), where
G has even-parity phase one and odd-parity phase t=e^(-2ic). The four masks
join one wire in each terminal pair. Nonvanishing on |t|=1 implies full
Schmidt rank for every real angle c on that cut.
"""

import json
from pathlib import Path

import sympy as s


ROOT = Path(__file__).resolve().parent
MASKS = (5, 6, 9, 10)
CUTS = ((0, 2), (0, 3))


def bit(value, wire):
    return (value >> (3 - wire)) & 1


def realign(unitary, cut):
    other = tuple(q for q in range(4) if q not in cut)
    result = s.zeros(16)
    for output in range(16):
        for input_ in range(16):
            oa = 2 * bit(output, cut[0]) + bit(output, cut[1])
            ia = 2 * bit(input_, cut[0]) + bit(input_, cut[1])
            ob = 2 * bit(output, other[0]) + bit(output, other[1])
            ib = 2 * bit(input_, other[0]) + bit(input_, other[1])
            result[4 * oa + ia, 4 * ob + ib] = unitary[output, input_]
    return result


def main():
    a = s.sqrt(2) / 2
    f = s.Matrix([[1, 0, 0, 0], [0, a, s.I * a, 0],
                  [0, s.I * a, a, 0], [0, 0, 0, 1]])
    suffix = s.kronecker_product(f, f)
    t = s.symbols("t")
    rows = []
    for mask in MASKS:
        gauge = s.diag(*[t if (basis & mask).bit_count() % 2 else 1
                         for basis in range(16)])
        u = gauge * suffix
        for cut in CUTS:
            factor = s.factor(realign(u, cut).det(method="domain-ge"))
            row = {"mask": mask, "cut": list(cut), "factor": str(factor)}
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    path = ROOT / "algebraic13_crosspair_two_body_symbolic_result.json"
    path.write_text(json.dumps({"row_phase_definition":
                                "even parity 1, odd parity t=e^(-2ic)",
                                "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
