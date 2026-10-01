"""Scan a two-qubit eigenspace gauge after output-CNOT relabelings.

For twelve ordered two-CNOT relabelings, the exact V4 eigenvalue labels
of local |01> and |10> on wires (2,3) agree for every spectator setting
on wires (0,1). Consequently RXX(theta) RYY(theta) on (2,3) changes the
eigenbasis but preserves diagonalization. Search exact pi-fraction theta
values and jointly resynthesize the final five-CNOT tail.
"""

import argparse
import json
import math
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from exact_check import exact_eigenvalue_labels
from qiskit_crosscheck import from_qiskit, to_qiskit
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent
ANGLES = (math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi)


def cx_index(x, edge):
    c, t = edge
    return x ^ (1 << (3 - t)) if x & (1 << (3 - c)) else x


def relabeled_labels(labels, first, second):
    # P = CX(second) CX(first), so P^-1 = CX(first) CX(second).
    return [labels[cx_index(cx_index(y, second), first)] for y in range(16)]


def gauge_condition(labels):
    # Spectator bits q0,q1; compare local q2q3=01 (index 1) and 10 (index 2).
    return all(labels[(s << 2) | 1] == labels[(s << 2) | 2]
               for s in range(4))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    split = 13 + len(exact_block()) + len(reverse_block())
    prefix, tail = source["gates"][:split], source["gates"][split:]
    edges = [(c, t) for c in range(4) for t in range(4) if c != t]
    compatible = [(first, second) for first in edges for second in edges
                  if first != second and gauge_condition(
                      relabeled_labels(labels, first, second))]
    assert len(compatible) == 12, compatible
    records = []
    best = None
    best_candidate = None
    for first, second in compatible:
        for angle in ANGLES:
            extras = [{"gate": "cx", "control": c, "target": t}
                      for c, t in (first, second)]
            circuit = to_qiskit({"n": 4, "gates": tail + extras})
            # Physical wires (2,3) are Qiskit qubits (1,0).
            circuit.rxx(angle, 1, 0)
            circuit.ryy(angle, 1, 0)
            optimized = transpile(circuit, basis_gates=["u", "cx"],
                                  optimization_level=3,
                                  seed_transpiler=args.seed)
            candidate = {"n": 4, "gates": prefix + from_qiskit(optimized)["gates"]}
            check = evaluate(candidate)
            record = {"output_cx_pair": [first, second],
                      "givens_angle_over_pi": round(angle / math.pi, 6),
                      "tail_cnot_count": optimized.count_ops().get("cx", 0),
                      "total_cnot_count": check["cnot_count"],
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid": check["valid_diagonalizer"]}
            records.append(record)
            if best is None or (record["total_cnot_count"],
                                record["off_diagonal_error"]) < (
                                    best["total_cnot_count"],
                                    best["off_diagonal_error"]):
                best, best_candidate = record, candidate
    args.output.write_text(json.dumps(best_candidate, indent=2) + "\n")
    counts = {}
    for record in records:
        count = record["tail_cnot_count"]
        counts[count] = counts.get(count, 0) + 1
    print(json.dumps({"compatible_output_pairs": compatible,
                      "label_condition": "D'(s,01)=D'(s,10) for all s in 00,01,10,11",
                      "trials": len(records), "tail_cnot_histogram": counts,
                      "best": best, "candidate": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
