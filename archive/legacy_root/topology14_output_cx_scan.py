"""Scan one output CNOT appended to the exact 15-CNOT Fourier tail.

An output CNOT permutes diagonal entries, so any successful transpiled
14-CNOT candidate remains a V4 diagonalizer. This bounded scan tries the
12 directed output CNOTs and retains the best resulting full gate list.
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    split = 13 + len(exact_block()) + len(reverse_block())
    prefix, tail = source["gates"][:split], source["gates"][split:]
    records = []
    best = None
    best_candidate = None
    for control in range(4):
        for target in range(4):
            if control == target:
                continue
            output_cx = {"gate": "cx", "control": control, "target": target}
            inp = {"n": 4, "gates": [*tail, output_cx]}
            optimized = transpile(to_qiskit(inp), basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=args.seed)
            emitted = from_qiskit(optimized)
            candidate = {"n": 4, "gates": prefix + emitted["gates"]}
            check = evaluate(candidate)
            record = {"output_cx": [control, target],
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
    print(json.dumps({"results": records, "best": best,
                      "candidate": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
