"""Jointly synthesize a full-parity-gauged Fourier boundary block.

R_P(delta), P=Z0Z1Z2Z3, commutes with V4. For the exact 15-CNOT open parity
prefix endpoint C=CNOT(3,1) CNOT(0,2), C P C=Z1Z2. Therefore inserting
R_{Z1Z2}(delta) before the inverse endpoint CNOTs corresponds exactly to
applying a commuting full-parity gauge on the input. Jointly transpile the
two inverse endpoint CNOTs, the two first Fourier butterflies, and this
cross-pair rotation. Grid delta=k*pi/16, k=0..16.
"""

import json
import math
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent


def boundary_block(baseline, delta):
    gates = []
    if delta:
        gates += [
            {"gate": "cx", "control": 1, "target": 2},
            {"gate": "rz", "qubit": 2, "theta": delta},
            {"gate": "cx", "control": 1, "target": 2},
        ]
    gates += [
        {"gate": "cx", "control": 0, "target": 2},
        {"gate": "cx", "control": 3, "target": 1},
    ]
    gates += baseline["gates"][18:38]
    return gates


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    exact15 = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix = exact15["gates"][:13]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    tail = baseline["gates"][38:]
    results = []
    best = None
    best_candidate = None
    for k in range(17):
        delta = k * math.pi / 16
        raw_block = {"n": 4, "gates": boundary_block(baseline, delta)}
        optimized = transpile(to_qiskit(raw_block),
                              basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        candidate = {"n": 4, "gates": prefix
                     + from_qiskit(optimized)["gates"] + tail}
        check = evaluate(candidate)
        record = {"k": k, "delta_pi_over_16": k,
                  "raw_boundary_cnot": sum(g["gate"] == "cx"
                                           for g in raw_block["gates"]),
                  "optimized_boundary_cnot": optimized.count_ops().get("cx", 0),
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
        path = ROOT / "algebraic14_full_parity_boundary_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        best["candidate"] = path.name
    print(json.dumps({"results": results, "best": best}, indent=2))


if __name__ == "__main__":
    main()
