"""Exact 15-CNOT V4 candidate from a six-CNOT parity prefix.

The open prefix has inverse CNOTs (3,1),(0,2). Each inverse CNOT fuses
with its disjoint two-qubit Fourier butterfly, producing a two-CNOT block.
Both replacement blocks use the exact Weyl entangler
RXX(-pi/2) RYY(-pi/4), implemented with two CNOTs and pi-fraction locals.
"""

import json
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix
from topology16_exact_block import exact_block


ROOT = Path(__file__).resolve().parent


def rot(name, q, theta):
    return {"gate": name, "qubit": q, "theta": theta}


def reverse_block():
    # Qiskit Weyl factors for CX(3,1) followed by baseline F2(1,3).
    # Qiskit local q0 is physical 3; local q1 is physical 1.
    c = {"gate": "cx", "control": 1, "target": 3}
    return [
        rot("rz", 3, "-pi"), rot("rz", 1, "-pi/2"),
        rot("ry", 3, "3*pi/4"), rot("ry", 1, "pi/2"),
        rot("rz", 1, "pi/2"),
        rot("rx", 3, "-pi/2"), rot("rx", 1, "-pi/2"),
        c,
        rot("rx", 1, "-pi/2"), rot("rz", 3, "-pi/4"),
        c,
        rot("rx", 3, "pi/2"), rot("rx", 1, "pi/2"),
        rot("rz", 3, "-pi"), rot("rz", 1, "3*pi/4"),
        rot("ry", 3, "pi/2"), rot("ry", 1, "pi"),
        rot("rz", 3, "-pi/2"), rot("rz", 1, "-pi/4"),
    ]


def local_block_unitary(gates, upper, lower):
    result = np.eye(4, dtype=complex)
    for g in gates:
        if g["gate"] == "cx":
            matrix = cnot(0 if g["control"] == upper else 1,
                          0 if g["target"] == upper else 1, 2)
        else:
            q = 0 if g["qubit"] == upper else 1
            matrix = embedded_one_qubit(one_qubit_matrix(g), q, 2)
        result = matrix @ result
    return result


def error_up_to_phase(left, right):
    phase = np.vdot(left.ravel(), right.ravel())
    return float(np.max(abs(left * phase / abs(phase) - right)))


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    numeric = json.loads((ROOT / "algebraic16_to15_gauge_1_best.json").read_text())
    prefix = numeric["gates"][:13]
    source_first = [{"gate": "cx", "control": 0, "target": 2},
                    *baseline["gates"][18:28]]
    source_second = [{"gate": "cx", "control": 3, "target": 1},
                     *baseline["gates"][28:38]]
    first, second = exact_block(), reverse_block()
    first_error = error_up_to_phase(local_block_unitary(first, 0, 2),
                                    local_block_unitary(source_first, 0, 2))
    second_error = error_up_to_phase(local_block_unitary(second, 1, 3),
                                     local_block_unitary(source_second, 1, 3))
    candidate = {"n": 4, "gates": prefix + first + second + baseline["gates"][38:]}
    output = ROOT / "topology16_15_exact_candidate.json"
    output.write_text(json.dumps(candidate, indent=2) + "\n")
    print(json.dumps({"first_block_error_up_to_phase": first_error,
                      "second_block_error_up_to_phase": second_error,
                      "cnot_count": sum(g["gate"] == "cx" for g in candidate["gates"]),
                      "output": str(output)}))


if __name__ == "__main__":
    main()
