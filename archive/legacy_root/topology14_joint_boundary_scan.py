"""Resynthesize the reverse butterfly jointly with the final Fourier tail.

The block acts on all four wires and contains seven CNOTs. Scan the identity
and each output CNOT, which can permute diagonal labels without changing the
diagonalization condition. A six-CNOT joint block yields 14 total.
"""

import argparse
import json
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit
from topology16_exact_block import exact_block


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    split = 13 + len(exact_block())
    prefix, joint = source["gates"][:split], source["gates"][split:]
    options = [None] + [(a, b) for a in range(4) for b in range(4) if a != b]
    records = []
    best = None
    best_candidate = None
    for option in options:
        extras = [] if option is None else [
            {"gate": "cx", "control": option[0], "target": option[1]}]
        inp = {"n": 4, "gates": joint + extras}
        optimized = transpile(to_qiskit(inp), basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=args.seed)
        emitted = from_qiskit(optimized)
        candidate = {"n": 4, "gates": prefix + emitted["gates"]}
        result = evaluate(candidate)
        record = {"output_cx": option,
                  "joint_cnot_count": optimized.count_ops().get("cx", 0),
                  "total_cnot_count": result["cnot_count"],
                  "off_diagonal_error": result["off_diagonal_error"],
                  "valid": result["valid_diagonalizer"]}
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
