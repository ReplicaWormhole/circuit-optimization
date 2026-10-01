"""Construct the rational-angle four-matchgate 14-CX diagonalizer.

The six-CX prefix is the verified exact prefix of the 15-CX circuit.  Four
two-CX Cartan blocks and Clifford/pi-six local frames implement its new tail.
Output diagonal Rz gates are omitted because they cannot affect diagonalization.
"""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from qiskit import transpile
from qiskit.quantum_info import Operator

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
PI = math.pi


def add_zyz(circuit, wire, theta, phi, lam):
    qubit = 3 - wire
    circuit.rz(lam * PI, qubit)
    circuit.ry(theta * PI, qubit)
    circuit.rz(phi * PI, qubit)


def make_circuit():
    prefix = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    circuit = to_qiskit({"n": 4, "gates": prefix["gates"][:13]})
    input_zyz = {
        0: (3 / 4, 0, 0),
        1: (1 / 2, -1 / 2, -1 / 3),
        2: (1 / 2, -3 / 2, 1 / 6),
        3: (3 / 4, -1, 1),
    }
    center_zyz = {
        0: (1 / 2, 0, -1),
        1: (1, 0, 3 / 2),
        2: (1, 0, 3 / 2),
        3: (1 / 2, -3 / 2, 0),
    }
    for wire in range(4):
        add_zyz(circuit, wire, *input_zyz[wire])
    for first, second in ((0, 2), (1, 3)):
        q1, q2 = 3 - first, 3 - second
        circuit.rxx(-PI / 2, q1, q2)
        circuit.ryy(-PI / 4, q1, q2)
    for wire in range(4):
        add_zyz(circuit, wire, *center_zyz[wire])
    for first, second in ((0, 1), (2, 3)):
        q1, q2 = 3 - first, 3 - second
        circuit.rxx(-PI / 4, q1, q2)
        circuit.ryy(-PI / 4, q1, q2)
    return circuit


def angle_string(value):
    value = Fraction(value)
    if value == 0:
        return "0*pi"
    sign = "-" if value < 0 else ""
    numerator, denominator = abs(value.numerator), value.denominator
    num = "" if numerator == 1 else str(numerator) + "*"
    den = "" if denominator == 1 else "/" + str(denominator)
    return f"{sign}{num}pi{den}"


def append_rotation(gates, name, wire, value):
    if value:
        gates.append({"gate": name, "qubit": wire,
                      "theta": angle_string(value)})


def append_zyz(gates, wire, theta, phi, lam):
    append_rotation(gates, "rz", wire, lam)
    append_rotation(gates, "ry", wire, theta)
    append_rotation(gates, "rz", wire, phi)


def append_cartan_two_cx(gates, first, second, a, b):
    """Implement exp(i*pi*(a XX+b YY)) exactly with two CNOTs.

    B=Rx(pi/2) tensor Rx(pi/2) fixes XX and maps ZZ to YY.  Conjugating
    Rx_first(-2a*pi) and Rz_second(-2b*pi) by CX gives XX and ZZ.
    """
    half = Fraction(1, 2)
    for wire in (first, second):
        append_rotation(gates, "rx", wire, -half)
    gates.append({"gate": "cx", "control": first, "target": second})
    append_rotation(gates, "rx", first, -2 * a)
    append_rotation(gates, "rz", second, -2 * b)
    gates.append({"gate": "cx", "control": first, "target": second})
    for wire in (first, second):
        append_rotation(gates, "rx", wire, half)


def make_rational_candidate():
    prefix = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    gates = list(prefix["gates"][:13])
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
    for wire in range(4):
        append_zyz(gates, wire, *input_zyz[wire])
    for first, second in ((0, 2), (1, 3)):
        append_cartan_two_cx(gates, first, second, Fraction(1, 4), Fraction(1, 8))
    for wire in range(4):
        append_zyz(gates, wire, *center_zyz[wire])
    for first, second in ((0, 1), (2, 3)):
        append_cartan_two_cx(gates, first, second, Fraction(1, 8), Fraction(1, 8))
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--rational-candidate", type=Path)
    args = parser.parse_args()
    circuit = make_circuit()
    raw = Operator(circuit).data
    optimized = transpile(circuit, basis_gates=["u", "cx"],
                          optimization_level=3, seed_transpiler=42)
    candidate = from_qiskit(optimized)
    args.candidate.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    qmat = circuit_matrix(candidate["gates"])
    overlap = (raw.conj() * qmat).sum()
    phase = overlap / abs(overlap)
    matrix_error = float(abs(raw - qmat / phase).max())
    result = {"candidate": str(args.candidate),
              "raw_circuit_cnot_count": circuit.count_ops().get("cx", 0),
              "transpiled_cnot_count": check["cnot_count"],
              "off_diagonal_error": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"],
              "raw_to_transpiled_matrix_error": matrix_error,
              "proof_status": "numerical pending exact symbolic verification"}
    if args.rational_candidate:
        rational = make_rational_candidate()
        args.rational_candidate.write_text(json.dumps(rational, indent=2) + "\n")
        rat_check = evaluate(rational)
        rat_matrix = circuit_matrix(rational["gates"])
        rat_overlap = (raw.conj() * rat_matrix).sum()
        rat_phase = rat_overlap / abs(rat_overlap)
        result["rational_candidate"] = str(args.rational_candidate)
        result["rational_cnot_count"] = rat_check["cnot_count"]
        result["rational_off_diagonal_error"] = rat_check["off_diagonal_error"]
        result["raw_to_rational_matrix_error"] = float(abs(raw - rat_matrix / rat_phase).max())
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
