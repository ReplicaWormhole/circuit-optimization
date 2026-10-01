"""Search exact CNOT-identity insertions at Fourier block boundaries.

At a boundary insert E;E for each of twelve directed CNOTs. Synthesize the
left block plus E and E plus the right block independently. The product is
exactly the original seven-CNOT joint block, up to transpiler precision.
"""

import argparse
import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent


def synth(gates, seed):
    circuit = transpile(to_qiskit({"n": 4, "gates": gates}),
                        basis_gates=["u", "cx"],
                        optimization_level=3, seed_transpiler=seed)
    return from_qiskit(circuit)["gates"], circuit.count_ops().get("cx", 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    fixed = 13 + len(exact_block())
    prefix, joint = source["gates"][:fixed], source["gates"][fixed:]
    assert len(joint) == len(reverse_block()) + 24
    splits = {"after_reverse_butterfly": len(reverse_block()),
              "after_middle_cz": len(reverse_block()) + 3,
              "after_first_final_f2": len(reverse_block()) + 15}
    records = []
    best = None
    best_candidate = None
    for name, split in splits.items():
        left, right = joint[:split], joint[split:]
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                edge = {"gate": "cx", "control": control, "target": target}
                left_gates, left_cx = synth(left + [edge], args.seed)
                right_gates, right_cx = synth([edge] + right, args.seed)
                candidate = {"n": 4, "gates": prefix + left_gates + right_gates}
                check = evaluate(candidate)
                record = {"boundary": name, "inserted_cx": [control, target],
                          "left_cnot_count": left_cx,
                          "right_cnot_count": right_cx,
                          "joint_cnot_count": left_cx + right_cx,
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
        count = record["joint_cnot_count"]
        counts[count] = counts.get(count, 0) + 1
    print(json.dumps({"trials": len(records), "joint_cnot_histogram": counts,
                      "best": best, "candidate": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
