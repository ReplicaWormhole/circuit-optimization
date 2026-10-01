"""Exact first-pair absorption screen for five-CNOT parity prefixes.

Enumerate every ordered five-CNOT circuit (12**5). A parity phase model with
five nonlocal masks requires that all five distinct masks appear as a target
row. Compare the final linear map with the exact six-CNOT prefix after one
CNOT on either first matchgate pair. The second check uses Q(zeta_96) to
classify the hypothetical first block with its exact input local frame.
"""

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parent
BASIS = (8, 4, 2, 1)
EXACT6_ROWS = (8, 5, 10, 1)
PAIRS = ((0, 2), (1, 3))
LOCAL = {
    0: (Fraction(3, 4), Fraction(0), Fraction(0)),
    1: (Fraction(1, 2), Fraction(-1, 2), Fraction(-1, 3)),
    2: (Fraction(1, 2), Fraction(-3, 2), Fraction(1, 6)),
    3: (Fraction(3, 4), Fraction(-1), Fraction(1)),
}


def endpoint_screen():
    models = json.loads((ROOT / "algebraic14_fivemask_models.json").read_text())
    supports = {frozenset(record["support"]) for record in models}
    assert len(models) == len(supports) == 54
    targets = {}
    for pair in PAIRS:
        for control, target in (pair, pair[::-1]):
            rows = list(EXACT6_ROWS)
            rows[target] ^= rows[control]
            targets[tuple(rows)] = f"{control}->{target}"

    endpoint_counts = Counter()
    catalog_counts = Counter()
    valid_paths = 0
    valid_endpoints = set()

    def walk(rows, masks, depth):
        nonlocal valid_paths
        if depth == 5:
            key = targets.get(rows)
            if key:
                endpoint_counts[key] += 1
            if masks in supports:
                valid_paths += 1
                valid_endpoints.add(rows)
                if key:
                    catalog_counts[key] += 1
            return
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                next_rows = list(rows)
                next_rows[target] ^= rows[control]
                mask = next_rows[target]
                next_masks = masks
                if mask.bit_count() > 1:
                    next_masks = masks | {mask}
                walk(tuple(next_rows), next_masks, depth + 1)

    walk(BASIS, frozenset(), 0)
    return {
        "model_count": len(models),
        "ordered_five_cnot_schedules": 12**5,
        "valid_model_paths": valid_paths,
        "valid_model_endpoints": len(valid_endpoints),
        "one_cnot_targets": {v: list(k) for k, v in targets.items()},
        "all_path_target_counts": dict(endpoint_counts),
        "valid_model_target_counts": {key: catalog_counts[key] for key in targets.values()},
    }


def cartan_screen():
    field = sp.QQ.cyclotomic_field(96)
    z = field.unit
    zero, one = field.zero, field.one
    half = field.from_sympy(sp.Rational(1, 2))
    imag = z**24
    h = (z**12 + z**-12) * half

    def matmul(a, b):
        return [[sum((a[i][k] * b[k][j] for k in range(len(b))), zero)
                 for j in range(len(b[0]))] for i in range(len(a))]

    def transpose(a):
        return [list(row) for row in zip(*a)]

    def kron(a, b):
        return [[a[i // 2][j // 2] * b[i % 2][j % 2] for j in range(4)]
                for i in range(4)]

    def identity(n):
        return [[one if i == j else zero for j in range(n)] for i in range(n)]

    def addscaled(a, b, scale):
        return [[a[i][j] + scale * b[i][j] for j in range(len(a[0]))]
                for i in range(len(a))]

    def trig(theta):
        exp = 48 * theta
        assert exp.denominator == 1
        positive, negative = z**int(exp), z**-int(exp)
        return (positive + negative) * half, (positive - negative) * half / imag

    def rz(theta):
        exp = 24 * theta
        assert exp.denominator == 1
        return [[z**-int(exp), zero], [zero, z**int(exp)]]

    def ry(theta):
        exp = 24 * theta
        assert exp.denominator == 1
        positive, negative = z**int(exp), z**-int(exp)
        co = (positive + negative) * half
        si = (positive - negative) * half / imag
        return [[co, -si], [si, co]]

    def zyz(angles):
        theta, phi, lam = angles
        return matmul(matmul(rz(phi), ry(theta)), rz(lam))

    x = [[zero, one], [one, zero]]
    y = [[zero, -imag], [imag, zero]]
    xx, yy = kron(x, x), kron(y, y)
    eye4 = identity(4)
    ca, sa = trig(Fraction(1, 4))
    cb, sb = trig(Fraction(1, 8))
    f = matmul(addscaled([[ca * v for v in row] for row in eye4], xx, imag * sa),
               addscaled([[cb * v for v in row] for row in eye4], yy, imag * sb))
    q = [[h, imag*h, zero, zero], [zero, zero, imag*h, h],
         [zero, zero, imag*h, -h], [h, -imag*h, zero, zero]]
    qdag = transpose([[h, -imag*h, zero, zero], [zero, zero, -imag*h, h],
                      [zero, zero, -imag*h, -h], [h, imag*h, zero, zero]])
    cx_forward = [[one, zero, zero, zero], [zero, one, zero, zero],
                  [zero, zero, zero, one], [zero, zero, one, zero]]
    cx_reverse = [[one, zero, zero, zero], [zero, zero, zero, one],
                  [zero, zero, one, zero], [zero, one, zero, zero]]

    rows = {}
    for pair in PAIRS:
        frame = kron(zyz(LOCAL[pair[0]]), zyz(LOCAL[pair[1]]))
        for name, cx in ((f"{pair[0]}->{pair[1]}", cx_forward),
                         (f"{pair[1]}->{pair[0]}", cx_reverse)):
            u = matmul(matmul(f, frame), cx)
            b = matmul(matmul(qdag, u), q)
            m = matmul(transpose(b), b)
            # det(U)=-1; multiply U by zeta_8=e^(i*pi/4), hence M by i,
            # to obtain SU(4)-normalized magic-basis invariants.
            normalized = [[imag * v for v in row] for row in m]
            traces = []
            power = identity(4)
            for _ in range(3):
                power = matmul(power, normalized)
                traces.append(sum((power[i][i] for i in range(4)), zero))
            rows[name] = {
                "trace_m_exact": [str(field.to_sympy(t)) for t in traces],
                "trace_m_zero": [t == zero for t in traces],
                "same_magic_polynomial_as_first_matchgate": all(t == zero for t in traces),
            }
    return rows


def main():
    result = {"endpoint_screen": endpoint_screen(), "cartan_screen": cartan_screen()}
    path = ROOT / "search13_firstpair_prefix_screen_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
