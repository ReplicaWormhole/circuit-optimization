"""Build an exact-angle 16-CNOT cycle diagonalizer.

The seven-CNOT parity prefix is from the endpoint (0,2) construction.
Replace the numeric KAK decomposition of CX(0,2) followed by the first
Fourier two-qubit block by an exact two-CNOT XX/YY decomposition:

RXX(-pi/2) RYY(-pi/4)
 = K CX(0,2) [RX_0(-pi/2) RZ_2(-pi/4)] CX(0,2) K†,
K = RX_0(pi/2) RX_2(pi/2).

The outer one-qubit rotations are exact Weyl factors for this block.
This formula is verified against the source block up to global phase.
"""

import json
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix


ROOT = Path(__file__).resolve().parent


def one(name, q, theta):
    return {"gate": name, "qubit": q, "theta": theta}


def exact_block():
    c = {"gate": "cx", "control": 0, "target": 2}
    return [
        one("rz", 2, "pi/2"), one("rz", 0, "pi"),
        one("ry", 2, "pi/2"), one("ry", 0, "pi/4"),
        one("rz", 2, "-pi/2"),
        one("rx", 2, "-pi/2"), one("rx", 0, "-pi/2"),
        c,
        one("rx", 0, "-pi/2"), one("rz", 2, "-pi/4"),
        c,
        one("rx", 2, "pi/2"), one("rx", 0, "pi/2"),
        one("rz", 2, "-pi"),
        one("rz", 0, "-pi"), one("ry", 0, "pi/2"),
        one("rz", 0, "pi/2"),
    ]


def two_qubit_unitary(gates):
    u = np.eye(4, dtype=complex)
    for gate in gates:
        if gate["gate"] == "cx":
            m = cnot(0, 1, 2)
        else:
            m = embedded_one_qubit(one_qubit_matrix(gate),
                                   0 if gate["qubit"] == 0 else 1, 2)
        u = m @ u
    return u


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    endpoint = json.loads((ROOT / "algebraic16_endpoint_0_2.json").read_text())
    prefix = endpoint["gates"][:14]
    original_block = [{"gate": "cx", "control": 0, "target": 2},
                      *baseline["gates"][18:28]]
    exact = exact_block()
    left, right = two_qubit_unitary(exact), two_qubit_unitary(original_block)
    phase = np.vdot(left.ravel(), right.ravel())
    error = np.max(np.abs(left * phase / abs(phase) - right))
    output = {"n": 4, "gates": prefix + exact + baseline["gates"][28:]}
    output_path = ROOT / "topology16_exact_angle_candidate.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"block_error_up_to_phase": float(error),
                      "cnot_count": sum(g["gate"] == "cx" for g in output["gates"]),
                      "output": str(output_path)}))


if __name__ == "__main__":
    main()
