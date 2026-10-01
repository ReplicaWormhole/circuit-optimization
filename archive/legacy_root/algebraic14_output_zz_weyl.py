"""Scan output ZZ gauges on the final two-qubit Fourier butterfly.

Left multiplication of a V4 diagonalizer by any computational-basis
diagonal unitary preserves diagonalization. For the final H_lift(2,3), test
RZZ(2,3,delta) H_lift(2,3), delta=k*pi/16, k=0..16. Report its Weyl
coordinates, Qiskit CNOT synthesis count, and the full-circuit residual.
"""

import json
import math
from pathlib import Path

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Operator
from qiskit.synthesis import TwoQubitWeylDecomposition

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent


def renumber(gates):
    result = []
    for gate in gates:
        item = gate.copy()
        if item["gate"] == "cx":
            item["control"] -= 2
            item["target"] -= 2
        else:
            item["qubit"] -= 2
        result.append(item)
    return result


def main():
    exact15 = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    butterfly = baseline["gates"][52:62]
    assert exact15["gates"][-10:] == butterfly
    prefix = exact15["gates"][:-10]
    results = []
    best = None
    best_candidate = None
    for k in range(17):
        delta = k * math.pi / 16
        gates = renumber(butterfly)
        if delta:
            gates += [
                {"gate": "cx", "control": 0, "target": 1},
                {"gate": "rz", "qubit": 1, "theta": delta},
                {"gate": "cx", "control": 0, "target": 1},
            ]
        qc = to_qiskit({"n": 2, "gates": gates})
        unitary = Operator(qc).data
        schmidt = unitary.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
        schmidt_rank = int(np.linalg.matrix_rank(schmidt, tol=1e-9))
        weyl = TwoQubitWeylDecomposition(unitary)
        optimized = transpile(qc, basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        emitted = from_qiskit(optimized)
        # Two-qubit indices 0,1 map back to physical 2,3.
        mapped = []
        for gate in emitted["gates"]:
            item = gate.copy()
            if item["gate"] == "cx":
                item["control"] += 2
                item["target"] += 2
            else:
                item["qubit"] += 2
            mapped.append(item)
        candidate = {"n": 4, "gates": prefix + mapped}
        check = evaluate(candidate)
        record = {"k": k, "delta_pi_over_16": k,
                  "weyl_a": float(weyl.a),
                  "weyl_b": float(weyl.b),
                  "weyl_c": float(weyl.c),
                  "operator_schmidt_rank": schmidt_rank,
                  "block_cnot": optimized.count_ops().get("cx", 0),
                  "total_cnot": check["cnot_count"],
                  "off_diagonal_error": check["off_diagonal_error"],
                  "valid": check["valid_diagonalizer"]}
        results.append(record)
        if best is None or (record["total_cnot"],
                            record["off_diagonal_error"]) < (
                                best["total_cnot"],
                                best["off_diagonal_error"]):
            best, best_candidate = record, candidate
    if best_candidate:
        path = ROOT / "algebraic14_output_zz_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        best["candidate"] = path.name
    print(json.dumps({"results": results, "best": best}, indent=2))


if __name__ == "__main__":
    main()
