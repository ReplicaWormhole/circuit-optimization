"""Export and symbolically certify an exact 16-CNOT V4 diagonalizer.

The only angle that is not a rational multiple of pi is
phi = pi - atan(1/sqrt(2)). SymPy verifies that the explicit two-CNOT
replacement equals the original two-CNOT Fourier butterfly preceded by
CNOT(0,2), up to global phase. The prefix parity-mask argument is documented
in algebraic16_exact_derivation.md.
"""

import json
import math
from pathlib import Path

import sympy as sp

from check_circuit import evaluate


ROOT = Path(__file__).parent
PHI = math.pi - math.atan(1 / math.sqrt(2))


def u(qubit, theta, phi, lam):
    return {"gate": "u3", "qubit": qubit,
            "theta": theta, "phi": phi, "lam": lam}


PREFIX = [
    {"gate": "rz", "qubit": 1, "theta": "5*pi/8"},
    {"gate": "rz", "qubit": 2, "theta": "3*pi/4"},
    {"gate": "rz", "qubit": 3, "theta": "3*pi/8"},
    {"gate": "cx", "control": 0, "target": 3},
    {"gate": "cx", "control": 1, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "-pi/2"},
    {"gate": "cx", "control": 1, "target": 3},
    {"gate": "rz", "qubit": 3, "theta": "pi/4"},
    {"gate": "cx", "control": 2, "target": 3},
    {"gate": "rz", "qubit": 3, "theta": "pi/8"},
    {"gate": "cx", "control": 0, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "3*pi/8"},
    {"gate": "cx", "control": 1, "target": 2},
    {"gate": "cx", "control": 2, "target": 3},
]


# Chronological two-CNOT implementation of H_lift(0,2) CNOT(0,2).
REPLACEMENT = [
    u(2, "pi/2", 0, "-pi"),
    u(0, "pi/2", "-pi/2", "-pi"),
    {"gate": "cx", "control": 2, "target": 0},
    u(2, "3*pi/4", "-pi", "-pi/2"),
    u(0, "2*pi/3", PHI, -PHI),
    {"gate": "cx", "control": 2, "target": 0},
    u(2, "pi/2", 0, "-pi"),
    u(0, "3*pi/4", "pi/2", "pi/2"),
]


def numeric_candidate():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    assert baseline["gates"][18]["gate"] == "rz"
    assert baseline["gates"][28]["gate"] == "rz"
    return {"n": 4, "gates": PREFIX + REPLACEMENT + baseline["gates"][28:]}


def exact_angle(value):
    if isinstance(value, (int, float)):
        if value == PHI:
            return sp.pi - sp.atan(1 / sp.sqrt(2))
        if value == -PHI:
            return -sp.pi + sp.atan(1 / sp.sqrt(2))
        return sp.Integer(value)
    return sp.sympify(value.replace("pi", "pi"), locals={"pi": sp.pi})


def exact_u3(gate):
    theta, phi, lam = (exact_angle(gate[key]) for key in ("theta", "phi", "lam"))
    if gate["phi"] == PHI and gate["lam"] == -PHI:
        return sp.Matrix([[sp.Rational(1, 2), (sp.sqrt(2) + sp.I) / 2],
                          [(-sp.sqrt(2) + sp.I) / 2, sp.Rational(1, 2)]])
    c, s = sp.cos(theta / 2), sp.sin(theta / 2)
    return sp.Matrix([[c, -sp.exp(sp.I * lam) * s],
                      [sp.exp(sp.I * phi) * s, sp.exp(sp.I * (phi + lam)) * c]])


def exact_one(gate):
    if gate["gate"] == "u3":
        return exact_u3(gate)
    theta = exact_angle(gate["theta"])
    c, s = sp.cos(theta / 2), sp.sin(theta / 2)
    if gate["gate"] == "rx":
        return sp.Matrix([[c, -sp.I * s], [-sp.I * s, c]])
    if gate["gate"] == "rz":
        return sp.diag(sp.exp(-sp.I * theta / 2), sp.exp(sp.I * theta / 2))
    raise ValueError(gate["gate"])


def exact_two_wire(gates):
    """Matrix on physical wires (0,2), in that left-to-right order."""
    eye = sp.eye(2)
    cx01 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                      [0, 0, 0, 1], [0, 0, 1, 0]])
    cx10 = sp.Matrix([[1, 0, 0, 0], [0, 0, 0, 1],
                      [0, 0, 1, 0], [0, 1, 0, 0]])
    result = sp.eye(4)
    for gate in gates:
        if gate["gate"] == "cx":
            matrix = cx01 if gate["control"] == 0 else cx10
        else:
            one = exact_one(gate)
            matrix = sp.kronecker_product(one, eye) if gate["qubit"] == 0 else sp.kronecker_product(eye, one)
        result = matrix * result
    return result


def certify_replacement():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    original = [{"gate": "cx", "control": 0, "target": 2}] + baseline["gates"][18:28]
    target = exact_two_wire(original)
    candidate = exact_two_wire(REPLACEMENT)
    entry = next((i, j) for i in range(4) for j in range(4) if target[i, j] != 0)
    phase = sp.simplify(candidate[entry] / target[entry])
    residual = candidate - phase * target
    checks = [sp.simplify(sp.expand_complex(entry)) == 0 for entry in residual]
    return bool(all(checks)), sp.simplify(phase)


def main():
    candidate = numeric_candidate()
    path = ROOT / "algebraic16_exact_candidate.json"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    symbolic_valid, global_phase = certify_replacement()
    print(json.dumps({"candidate": path.name,
                      "numeric": {k: v for k, v in evaluate(candidate).items()
                                  if k != "diagonal_eigenvalues"},
                      "symbolic_two_qubit_identity": symbolic_valid,
                      "symbolic_global_phase": str(global_phase)}, indent=2))
    if not symbolic_valid:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
