"""Rank 14-CNOT delete-and-reroute mutations by raw cycle residual.

This is a topology seed scan with one-qubit gates left at their exact
15-CNOT values. Its scores are only heuristics for subsequent optimization.
"""

import argparse
import json
from pathlib import Path

from check_circuit import evaluate


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-candidate", type=Path,
                        default=ROOT / "topology16_15_exact_candidate.json")
    parser.add_argument("--fixed-prefix-cx", type=int, default=0)
    args = parser.parse_args()
    baseline = json.loads(args.base_candidate.read_text())
    gates = baseline["gates"]
    cx_positions = [i for i, gate in enumerate(gates) if gate["gate"] == "cx"]
    records = []
    for deleted_cx, deleted_position in enumerate(cx_positions):
        if deleted_cx < args.fixed_prefix_cx:
            continue
        shortened = gates[:deleted_position] + gates[deleted_position + 1:]
        remaining = [i for i, gate in enumerate(shortened) if gate["gate"] == "cx"]
        for replace_cx, position in enumerate(remaining):
            if replace_cx < args.fixed_prefix_cx:
                continue
            old = shortened[position]
            for control in range(4):
                for target in range(4):
                    if control == target or (control, target) in (
                            (old["control"], old["target"]),
                            (old["target"], old["control"])):
                        continue
                    mutated = list(shortened)
                    mutated[position] = {"gate": "cx", "control": control,
                                         "target": target}
                    residual = evaluate({"n": 4, "gates": mutated})["off_diagonal_error"]
                    records.append({"deleted_cx": deleted_cx,
                                    "replace_cx": replace_cx,
                                    "old_edge": [old["control"], old["target"]],
                                    "new_edge": [control, target],
                                    "raw_offdiag": residual})
    records.sort(key=lambda r: r["raw_offdiag"])
    print(json.dumps({"base_candidate": str(args.base_candidate),
                      "fixed_prefix_cx": args.fixed_prefix_cx,
                      "mutations_evaluated": len(records),
                      "top_25": records[:25]}, indent=2))


if __name__ == "__main__":
    main()
