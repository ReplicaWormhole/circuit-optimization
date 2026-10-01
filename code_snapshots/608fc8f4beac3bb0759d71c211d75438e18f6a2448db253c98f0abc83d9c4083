"""Independently simulate and transpile a candidate with optional Qiskit."""

import argparse
import json
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

from check_circuit import angle, evaluate, right_shift


def to_qiskit(candidate):
    n = candidate["n"]
    circuit = QuantumCircuit(n)
    for gate in candidate["gates"]:
        name = gate["gate"].lower()
        # Qiskit q[0] is the least-significant bit in Operator's matrix.
        if name == "cx":
            circuit.cx(n - 1 - gate["control"], n - 1 - gate["target"])
        elif name in ("rx", "ry", "rz"):
            getattr(circuit, name)(angle(gate["theta"]), n - 1 - gate["qubit"])
        elif name == "u3":
            circuit.u(*(angle(gate[key]) for key in ("theta", "phi", "lam")),
                      n - 1 - gate["qubit"])
        else:
            getattr(circuit, name)(n - 1 - gate["qubit"])
    return circuit


def from_qiskit(circuit):
    n = circuit.num_qubits
    gates = []
    for instruction in circuit.data:
        name = instruction.operation.name
        qubits = [n - 1 - circuit.find_bit(qubit).index for qubit in instruction.qubits]
        if name == "cx":
            gates.append({"gate": "cx", "control": qubits[0], "target": qubits[1]})
        elif name == "u":
            theta, phi, lam = (float(value) for value in instruction.operation.params)
            gates.append({"gate": "u3", "qubit": qubits[0],
                          "theta": theta, "phi": phi, "lam": lam})
        else:
            raise ValueError(f"unexpected transpiled gate: {name}")
    return {"n": n, "gates": gates}


def off_diagonal_error(circuit):
    unitary = Operator(circuit).data
    conjugated = unitary @ right_shift(circuit.num_qubits) @ unitary.conj().T
    return float(np.max(np.abs(conjugated - np.diag(np.diag(conjugated)))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path,
                        help="save the transpiled gate list in checker format")
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text())
    original = to_qiskit(candidate)
    optimized = transpile(original, basis_gates=["u", "cx"],
                          optimization_level=3, seed_transpiler=42)
    emitted = from_qiskit(optimized)
    result = {
        "qiskit_input_cnot_count": original.count_ops().get("cx", 0),
        "qiskit_input_off_diagonal_error": off_diagonal_error(original),
        "qiskit_transpiled_cnot_count": optimized.count_ops().get("cx", 0),
        "qiskit_transpiled_off_diagonal_error": off_diagonal_error(optimized),
        "independent_checker": evaluate(emitted),
    }
    print(json.dumps(result, indent=2))
    if args.output:
        args.output.write_text(json.dumps(emitted, indent=2) + "\n")


if __name__ == "__main__":
    main()
