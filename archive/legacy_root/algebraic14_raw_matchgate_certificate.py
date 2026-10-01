"""Exact Q(zeta_96) certificate for the rational-angle 14-CX matchgate circuit.

The matrix is assembled directly from the six-CX prefix and the four XX/YY
rotation blocks, independently of Qiskit's floating-point KAK decomposition.
"""

import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi


ROOT = Path(__file__).resolve().parent
N = 96
F = sp.QQ.cyclotomic_field(N)
Z = F.unit
ZERO, ONE = F.zero, F.one
HALF = F.from_sympy(sp.Rational(1, 2))
I = Z ** (N // 4)


def exp_i_pi(angle):
    r = rational_pi(angle)
    assert N % (4 * r.denominator) == 0
    return Z ** (r.numerator * (N // (4 * r.denominator)))


def trig_half(angle):
    p = exp_i_pi(angle)
    return (p + p ** -1) * HALF, (p - p ** -1) / (2 * I)


def phase(angle):
    r = rational_pi(angle)
    assert N % (2 * r.denominator) == 0
    return Z ** (r.numerator * (N // (2 * r.denominator)))


def pi_fraction(value):
    r = sp.Rational(value).limit_denominator(48)
    assert r == value
    return f"{r.numerator}*pi/{r.denominator}"


def one_qubit(U, wire, matrix):
    mask = 1 << (3 - wire)
    for a in range(16):
        if a & mask:
            continue
        b = a ^ mask
        ra, rb = U[a], U[b]
        U[a] = [matrix[0][0] * x + matrix[0][1] * y for x, y in zip(ra, rb)]
        U[b] = [matrix[1][0] * x + matrix[1][1] * y for x, y in zip(ra, rb)]


def xx_yy(U, first, second, axis, angle):
    c, s = trig_half(angle)
    mask1, mask2 = 1 << (3 - first), 1 << (3 - second)
    mask = mask1 | mask2
    for a in range(16):
        b = a ^ mask
        if a > b:
            continue
        sign = ONE if axis == "xx" else (ONE if bool(a & mask1) != bool(a & mask2) else -ONE)
        ra, rb = U[a], U[b]
        U[a] = [c * x - I * s * sign * y for x, y in zip(ra, rb)]
        U[b] = [c * y - I * s * sign * x for x, y in zip(ra, rb)]


def apply_prefix(U, gate):
    name = gate["gate"].lower()
    if name == "cx":
        control, target = gate["control"], gate["target"]
        cm, tm = 1 << (3 - control), 1 << (3 - target)
        old = U[:]
        for a in range(16):
            U[a ^ tm if a & cm else a] = old[a]
        return
    w = gate["qubit"]
    if name == "u3":
        c, s = trig_half(gate["theta"])
        p, l = phase(gate["phi"]), phase(gate["lam"])
        M = ((c, -l * s), (p * s, p * l * c))
    elif name == "rz":
        p = exp_i_pi(gate["theta"])
        M = ((p ** -1, ZERO), (ZERO, p))
    elif name == "ry":
        c, s = trig_half(gate["theta"])
        M = ((c, -s), (s, c))
    elif name == "rx":
        c, s = trig_half(gate["theta"])
        M = ((c, -I * s), (-I * s, c))
    elif name == "h":
        h = (Z ** (N // 8) + Z ** (-N // 8)) * HALF
        M = ((h, h), (h, -h))
    elif name == "x":
        M = ((ZERO, ONE), (ONE, ZERO))
    else:
        raise ValueError(gate)
    one_qubit(U, w, M)


def add_zyz(U, wire, theta, phi, lam):
    for name, value in (("rz", lam), ("ry", theta), ("rz", phi)):
        apply_prefix(U, {"gate": name, "qubit": wire, "theta": pi_fraction(value)})


def main():
    prefix = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())["gates"][:13]
    U = [[ONE if a == b else ZERO for b in range(16)] for a in range(16)]
    for gate in prefix:
        apply_prefix(U, gate)
    inp = {0: (sp.Rational(3, 4), 0, 0),
           1: (sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(-1, 3)),
           2: (sp.Rational(1, 2), sp.Rational(-3, 2), sp.Rational(1, 6)),
           3: (sp.Rational(3, 4), -1, 1)}
    mid = {0: (sp.Rational(1, 2), 0, -1),
           1: (1, 0, sp.Rational(3, 2)),
           2: (1, 0, sp.Rational(3, 2)),
           3: (sp.Rational(1, 2), sp.Rational(-3, 2), 0)}
    for w in range(4):
        add_zyz(U, w, *inp[w])
    for a, b in ((0, 2), (1, 3)):
        xx_yy(U, a, b, "xx", "-pi/2")
        xx_yy(U, a, b, "yy", "-pi/4")
    for w in range(4):
        add_zyz(U, w, *mid[w])
    for a, b in ((0, 1), (2, 3)):
        xx_yy(U, a, b, "xx", "-pi/4")
        xx_yy(U, a, b, "yy", "-pi/4")
    roots = (ONE, I, -ONE, -I)
    labels = []
    for row in U:
        valid = [j for j, z in enumerate(roots)
                 if all(row[((b & 1) << 3) | (b >> 1)] == z * row[b] for b in range(16))]
        labels.append(valid)
    gate_list_matches = {}
    cnot_counts = {}
    for path in ("topology14_exact_matchgate_rational.json", "root14_exact_candidate.json"):
        gates = json.loads((ROOT / path).read_text())["gates"]
        V = [[ONE if a == b else ZERO for b in range(16)] for a in range(16)]
        for gate in gates:
            apply_prefix(V, gate)
        ratio = next(U[a][b] / V[a][b] for a in range(16) for b in range(16)
                     if V[a][b] != ZERO)
        gate_list_matches[path] = all(U[a][b] == ratio * V[a][b]
                                      for a in range(16) for b in range(16))
        cnot_counts[path] = sum(g["gate"] == "cx" for g in gates)
    result = {"field": "Q(zeta_96)", "exact_diagonalizer": all(len(x) == 1 for x in labels),
              "labels": [x[0] if len(x) == 1 else x for x in labels],
              "prefix_cnot_count": sum(g["gate"] == "cx" for g in prefix),
              "nonlocal_blocks": 4,
              "gate_list_cnot_counts": cnot_counts,
              "raw_equals_gate_lists_up_to_global_phase": gate_list_matches}
    (ROOT / "algebraic14_raw_matchgate_certificate_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["exact_diagonalizer"] and all(gate_list_matches.values())
                     and all(count == 14 for count in cnot_counts.values()) else 1)


if __name__ == "__main__":
    main()
