"""Synthesize the orbit-eigenbasis of the Bell-label signed swap.

The two Bell labels are a,b in {0,1,2,3}. After computing d=a xor b,
V4 acts as |a,d> -> (-1)^[a xor d=3] |a xor d,d>. Each d-sector is
diagonalized directly in two-state orbits, allowing all degenerate sectors.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix
from qiskit_crosscheck import from_qiskit


def index(a, d):
    return ((a >> 1) << 3) | ((d >> 1) << 2) | ((a & 1) << 1) | (d & 1)


def prefix_unitary():
    h = one_qubit_matrix({"gate": "h"})
    h0 = embedded_one_qubit(h, 0, 4)
    h1 = embedded_one_qubit(h, 1, 4)
    u = np.eye(16, dtype=complex)
    for g in [cnot(0, 2, 4), h0, cnot(1, 3, 4), h1,
              cnot(0, 1, 4), cnot(2, 3, 4)]:
        u = g @ u
    return u


def block_diagonalizer():
    u = np.zeros((16, 16), dtype=complex)
    for d in range(4):
        if d == 0:
            for a in range(4):
                u[index(a, d), index(a, d)] = 1
            continue
        vectors = []
        for a in range(4):
            b = a ^ d
            if a >= b:
                continue
            sign_a = -1 if a == 3 else 1
            sign_b = -1 if b == 3 else 1
            lam = 1j if sign_a * sign_b == -1 else 1
            for eigenvalue in (lam, -lam):
                v = np.zeros(4, dtype=complex)
                v[a] = 1 / math.sqrt(2)
                v[b] = eigenvalue / sign_a / math.sqrt(2)
                vectors.append(v)
        for out_a, v in enumerate(vectors):
            for in_a, coeff in enumerate(v):
                u[index(out_a, d), index(in_a, d)] = coeff.conjugate()
    return u


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    u = block_diagonalizer() @ prefix_unitary()
    circuit = QuantumCircuit(4)
    circuit.unitary(Operator(u), [3, 2, 1, 0])
    optimized = transpile(circuit, basis_gates=["u", "cx"],
                          optimization_level=3, seed_transpiler=args.seed)
    candidate = from_qiskit(optimized)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    matrix = Operator(optimized).data
    print(json.dumps({"unitarity_error": float(np.max(abs(u @ u.conj().T - np.eye(16)))),
                      "synthesis_matrix_error_up_to_global_phase":
                      float(np.max(abs(matrix * np.vdot(matrix.ravel(), u.ravel()) /
                                       abs(np.vdot(matrix.ravel(), u.ravel())) - u))),
                      "cnot_count": int(optimized.count_ops().get("cx", 0)),
                      "gate_count": len(candidate["gates"]),
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
