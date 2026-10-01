"""Test a shift-invariant opposite-pair diagonal gauge at the Fourier boundary.

J(delta)=R_{Z0Z2}(delta) R_{Z1Z3}(delta) commutes with V4. The exact15
prefix endpoint is C=CNOT(3,1) CNOT(0,2), and C J C becomes local
Rz(2,delta) Rz(1,delta). Thus insert those local rotations immediately
before the first two lifted Fourier butterflies and test whether their Weyl
class admits a smaller CNOT decomposition. Grid: delta=k*pi/16, k=0..16.
"""

import json
import math
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent


def optimized_block(inverse_cnot, local_qubit, original, delta):
    input_gates = []
    if delta:
        input_gates.append({"gate": "rz", "qubit": local_qubit,
                            "theta": delta})
    input_gates.append({"gate": "cx", "control": inverse_cnot[0],
                        "target": inverse_cnot[1]})
    input_gates.extend(original)
    optimized = transpile(to_qiskit({"n": 4, "gates": input_gates}),
                          basis_gates=["u", "cx"],
                          optimization_level=3, seed_transpiler=42)
    return from_qiskit(optimized)["gates"], optimized.count_ops().get("cx", 0)


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    exact15 = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix = exact15["gates"][:13]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    results = []
    best = None
    best_candidate = None
    for k in range(17):
        delta = k * math.pi / 16
        first, first_cx = optimized_block((0, 2), 2,
                                          baseline["gates"][18:28], delta)
        second, second_cx = optimized_block((3, 1), 1,
                                            baseline["gates"][28:38], delta)
        candidate = {"n": 4, "gates": prefix + first + second
                     + baseline["gates"][38:]}
        check = evaluate(candidate)
        record = {"k": k, "delta_pi_over_16": k,
                  "first_butterfly_cnot": first_cx,
                  "second_butterfly_cnot": second_cx,
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
        path = ROOT / "algebraic14_opposite_pair_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        best["candidate"] = path.name
    print(json.dumps({"results": results, "best": best}, indent=2))


if __name__ == "__main__":
    main()
