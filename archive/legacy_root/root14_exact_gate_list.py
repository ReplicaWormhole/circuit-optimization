"""Emit an exact rational-pi 14-CNOT circuit from four XX/YY blocks.

Identity used for any a,b:
  RXX(a) RYY(b) = K CX(c,t) [RX_c(a) RZ_t(b)] CX(c,t) K†,
  K = RX_c(pi/2) RX_t(pi/2).
The prefix and all remaining local rotations have rational-pi angles.
The output-only diagonal gates are omitted because they commute with D.
"""

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "root14_exact_candidate.json"


def pi(value):
    value = Fraction(value)
    if value == 0:
        return "0*pi"
    return f"{value.numerator}*pi/{value.denominator}"


def rotation(gates, name, qubit, angle):
    if Fraction(angle) != 0:
        gates.append({"gate": name, "qubit": qubit, "theta": pi(angle)})


def zyz(gates, qubit, theta, phi, lam):
    rotation(gates, "rz", qubit, lam)
    rotation(gates, "ry", qubit, theta)
    rotation(gates, "rz", qubit, phi)


def xx_yy(gates, control, target, a, b):
    rotation(gates, "rx", control, Fraction(-1, 2))
    rotation(gates, "rx", target, Fraction(-1, 2))
    gates.append({"gate": "cx", "control": control, "target": target})
    rotation(gates, "rx", control, a)
    rotation(gates, "rz", target, b)
    gates.append({"gate": "cx", "control": control, "target": target})
    rotation(gates, "rx", control, Fraction(1, 2))
    rotation(gates, "rx", target, Fraction(1, 2))


def make_candidate():
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    gates = source["gates"][:13].copy()
    input_zyz = {
        0: (Fraction(3, 4), 0, 0),
        1: (Fraction(1, 2), Fraction(-1, 2), Fraction(-1, 3)),
        2: (Fraction(1, 2), Fraction(-3, 2), Fraction(1, 6)),
        3: (Fraction(3, 4), -1, 1),
    }
    center_zyz = {
        0: (Fraction(1, 2), 0, -1),
        1: (1, 0, Fraction(3, 2)),
        2: (1, 0, Fraction(3, 2)),
        3: (Fraction(1, 2), Fraction(-3, 2), 0),
    }
    for q, angles in input_zyz.items():
        zyz(gates, q, *angles)
    for control, target in ((0, 2), (1, 3)):
        xx_yy(gates, control, target, Fraction(-1, 2), Fraction(-1, 4))
    for q, angles in center_zyz.items():
        zyz(gates, q, *angles)
    for control, target in ((0, 1), (2, 3)):
        xx_yy(gates, control, target, Fraction(-1, 4), Fraction(-1, 4))
    return {"n": 4, "gates": gates}


if __name__ == "__main__":
    candidate = make_candidate()
    OUTPUT.write_text(json.dumps(candidate, indent=2) + "\n")
    print(json.dumps({"candidate": str(OUTPUT),
                      "gates": len(candidate["gates"]),
                      "cnot_count": sum(g["gate"] == "cx"
                                        for g in candidate["gates"])}, indent=2))
