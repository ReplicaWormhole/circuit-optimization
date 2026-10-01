"""Operator-Schmidt ranks of the parity-gauged four-wire boundary.

A CNOT across a one-wire bipartition has Schmidt rank 2. If a circuit has
only three CNOTs on four wires, its interaction graph is a tree or smaller,
so at least two wires participate in at most one CNOT. Those wires have
operator-Schmidt rank at most 2, regardless of local gates. A target with
rank 4 on all four wires therefore cannot have a three-CNOT realization.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix


ROOT = Path(__file__).resolve().parent


def unitary(gates):
    result = np.eye(16, dtype=complex)
    for gate in gates:
        if gate["gate"] == "cx":
            matrix = cnot(gate["control"], gate["target"], 4)
        else:
            matrix = embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4)
        result = matrix @ result
    return result


def one_wire_schmidt_values(matrix, qubit):
    tensor = matrix.reshape((2,) * 8)
    axes = [qubit, 4 + qubit]
    axes += [q for q in range(4) if q != qubit]
    axes += [4 + q for q in range(4) if q != qubit]
    return np.linalg.svd(tensor.transpose(axes).reshape(4, 64), compute_uv=False)


def boundary(baseline, delta):
    gates = []
    if delta:
        gates += [{"gate": "cx", "control": 1, "target": 2},
                  {"gate": "rz", "qubit": 2, "theta": delta},
                  {"gate": "cx", "control": 1, "target": 2}]
    gates += [{"gate": "cx", "control": 0, "target": 2},
              {"gate": "cx", "control": 3, "target": 1}]
    gates += baseline["gates"][18:38]
    return gates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ks", type=int, nargs="+", default=list(range(17)))
    args = parser.parse_args()
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    rows = []
    for k in args.ks:
        matrix = unitary(boundary(baseline, k * math.pi / 16))
        spectra = [one_wire_schmidt_values(matrix, q) for q in range(4)]
        ranks = [int(np.count_nonzero(s > 1e-10)) for s in spectra]
        rows.append({"delta_pi_over_16": k, "one_wire_ranks": ranks,
                     "singular_values": [[float(s) for s in values]
                                         for values in spectra]})
    print(json.dumps({"tolerance": 1e-10, "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
