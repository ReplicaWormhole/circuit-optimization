"""Create a machine-readable exact 16-CNOT gate list and verify U V4 = D U.

The exact specification uses one named algebraic one-qubit gate. Its full
matrix is given in the JSON, and the verifier below constructs it directly
in Q(zeta_32). A separate checker-compatible file substitutes its U3 angles
numerically, using phi = pi - atan(1/sqrt(2)).
"""

import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp

from check_circuit import evaluate
from exact_check import rational_pi


ROOT = Path(__file__).parent
SPECIAL = {
    "gate": "u3_algebraic", "qubit": 0,
    "theta": "2*pi/3", "phi": "pi-atan(1/sqrt(2))",
    "lam": "-pi+atan(1/sqrt(2))",
    "matrix": [["1/2", "(sqrt(2)+i)/2"],
               ["(-sqrt(2)+i)/2", "1/2"]],
}


def exact_spec():
    numeric = json.loads((ROOT / "algebraic16_exact_candidate.json").read_text())
    gates = []
    special_count = 0
    for gate in numeric["gates"]:
        if gate["gate"] == "u3" and isinstance(gate.get("phi"), float):
            assert gate["qubit"] == 0 and gate["theta"] == "2*pi/3"
            gates.append(SPECIAL)
            special_count += 1
        else:
            gates.append({key: ("0*pi" if key in ("theta", "phi", "lam")
                                 and value == 0 else value)
                          for key, value in gate.items()})
    assert special_count == 1
    return {"n": 4, "gates": gates}


def numeric_from_spec(spec):
    phi = math.pi - math.atan(1 / math.sqrt(2))
    gates = []
    for gate in spec["gates"]:
        if gate["gate"] == "u3_algebraic":
            assert gate == SPECIAL
            gates.append({"gate": "u3", "qubit": 0,
                          "theta": "2*pi/3", "phi": phi, "lam": -phi})
        else:
            gates.append(gate)
    return {"n": 4, "gates": gates}


def exact_cycle_check(spec):
    if spec["n"] != 4:
        raise ValueError("only four qubits supported")
    field = sp.QQ.cyclotomic_field(32)
    zeta = field.unit
    zero, one = field.zero, field.one
    half = field.from_sympy(sp.Rational(1, 2))
    imaginary = zeta ** 8
    sqrt2 = zeta ** 4 + zeta ** -4
    hadamard_value = sqrt2 * half

    def phase(value, half_angle=False):
        frac = rational_pi(value)
        denominator = (4 if half_angle else 2) * frac.denominator
        if 32 % denominator:
            raise ValueError(f"angle {value} outside cyclotomic field")
        return zeta ** (frac.numerator * (32 // denominator))

    def local(gate):
        name = gate["gate"]
        if name == "u3_algebraic":
            if gate != SPECIAL:
                raise ValueError("unrecognized algebraic gate")
            return ((half, (sqrt2 + imaginary) * half),
                    ((-sqrt2 + imaginary) * half, half))
        if name == "h":
            h = hadamard_value
            return ((h, h), (h, -h))
        if name == "x":
            return ((zero, one), (one, zero))
        if name == "s":
            return ((one, zero), (zero, imaginary))
        if name == "sdg":
            return ((one, zero), (zero, -imaginary))
        if name in ("rx", "ry", "rz", "u3"):
            p = phase(gate["theta"], half_angle=True)
            q = p ** -1
            cosine = (p + q) * half
            sine = (p - q) / (2 * imaginary)
            if name == "rx":
                return ((cosine, -imaginary * sine),
                        (-imaginary * sine, cosine))
            if name == "ry":
                return ((cosine, -sine), (sine, cosine))
            if name == "rz":
                return ((q, zero), (zero, p))
            phi, lam = phase(gate["phi"]), phase(gate["lam"])
            return ((cosine, -lam * sine),
                    (phi * sine, phi * lam * cosine))
        raise ValueError(name)

    rows = [[one if i == j else zero for j in range(16)] for i in range(16)]
    cx_count = 0
    for position, gate in enumerate(spec["gates"]):
        if gate["gate"] == "cx":
            c, t = gate["control"], gate["target"]
            if not 0 <= c < 4 or not 0 <= t < 4 or c == t:
                raise ValueError(f"invalid CNOT at {position}")
            c_mask, t_mask = 1 << (3 - c), 1 << (3 - t)
            next_rows = [None] * 16
            for row in range(16):
                out = row ^ t_mask if row & c_mask else row
                next_rows[out] = rows[row]
            rows = next_rows
            cx_count += 1
        else:
            qubit = gate["qubit"]
            matrix = local(gate)
            mask = 1 << (3 - qubit)
            for row in range(16):
                if row & mask:
                    continue
                partner = row ^ mask
                first, second = rows[row], rows[partner]
                rows[row] = [matrix[0][0] * a + matrix[0][1] * b
                             for a, b in zip(first, second)]
                rows[partner] = [matrix[1][0] * a + matrix[1][1] * b
                                 for a, b in zip(first, second)]

    roots = (one, imaginary, -one, -imaginary)
    labels = []
    for row in rows:
        matches = [index for index, root in enumerate(roots)
                   if all(row[((basis & 1) << 3) | (basis >> 1)] == root * row[basis]
                          for basis in range(16))]
        if len(matches) != 1:
            return {"exact_diagonalizer": False, "cnot_count": cx_count,
                    "output_labels": labels}
        labels.append(matches[0])
    return {"exact_diagonalizer": True, "cnot_count": cx_count,
            "cyclotomic_conductor": 32,
            "output_labels": labels,
            "eigenvalue_multiplicities": [labels.count(k) for k in range(4)]}


def main():
    spec = exact_spec()
    spec_path = ROOT / "algebraic16_exact_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")
    numeric = numeric_from_spec(spec)
    numeric_path = ROOT / "algebraic16_from_exact_spec.json"
    numeric_path.write_text(json.dumps(numeric, indent=2) + "\n")
    check = exact_cycle_check(spec)
    result = {"exact_spec": spec_path.name,
              "numeric_export": numeric_path.name,
              "exact_check": check,
              "numeric_check": {k: v for k, v in evaluate(numeric).items()
                                if k != "diagonal_eigenvalues"}}
    print(json.dumps(result, indent=2))
    if not check["exact_diagonalizer"] or check["cnot_count"] != 16:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
