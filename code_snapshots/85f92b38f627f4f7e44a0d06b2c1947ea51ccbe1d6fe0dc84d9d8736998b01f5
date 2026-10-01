"""Synthesize explicit Bell-label two-qubit eigenbasis multiplexers.

For a=(a0,a1), b=(b0,b1), and d=a xor b, the conjugated cycle acts in
each d-sector as CZ_a X_a^d. The exact 4x4 eigenbasis is synthesized only
on the a register and controlled by the two bits of d.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Operator

from qiskit_crosscheck import from_qiskit


def sector_basis(d):
    if d == 0:
        return np.eye(4, dtype=complex)
    vectors = []
    for a in range(4):
        b = a ^ d
        if a >= b:
            continue
        s_a = -1 if a == 3 else 1
        s_b = -1 if b == 3 else 1
        lam = 1j if s_a * s_b == -1 else 1
        for z in (lam, -lam):
            v = np.zeros(4, dtype=complex)
            v[a] = 1 / math.sqrt(2)
            v[b] = z / s_a / math.sqrt(2)
            vectors.append(v)
    return np.array(vectors).conj()


def make_circuit():
    # Qiskit q[3],q[2],q[1],q[0] represent physical q0,q1,q2,q3.
    qc = QuantumCircuit(4)
    qc.cx(3, 1)
    qc.h(3)
    qc.cx(2, 0)
    qc.h(2)
    qc.cx(3, 2)
    qc.cx(1, 0)
    for d in range(1, 4):
        d0, d1 = d >> 1, d & 1
        if not d0:
            qc.x(2)
        if not d1:
            qc.x(0)
        gate = UnitaryGate(sector_basis(d)).control(2)
        # Controls: d1 q[0], d0 q[2]. Targets: a1 q[1], a0 q[3].
        qc.append(gate, [0, 2, 1, 3])
        if not d0:
            qc.x(2)
        if not d1:
            qc.x(0)
    return qc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    circuit = make_circuit()
    optimized = transpile(circuit, basis_gates=["u", "cx"],
                          optimization_level=3, seed_transpiler=args.seed)
    candidate = from_qiskit(optimized)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    print(json.dumps({"cnot_count": int(optimized.count_ops().get("cx", 0)),
                      "gate_count": len(candidate["gates"]),
                      "unitarity_error": float(np.max(abs(Operator(optimized).data @
                                         Operator(optimized).data.conj().T-np.eye(16)))),
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
