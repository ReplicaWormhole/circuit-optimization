"""Insert one directed CNOT into the best 15-CNOT Bell replacement."""

import argparse
import json
from pathlib import Path

from check_circuit import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidate = json.loads((Path(__file__).parent / "last_f2_swap_basis_best.json").read_text())
    gates = candidate["gates"]
    assert sum(gate["gate"] == "cx" for gate in gates) == 15
    best = None
    count = 0
    for position in range(len(gates) + 1):
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                patch = {"gate": "cx", "control": control, "target": target}
                trial = {"n": 4, "gates": gates[:position] + [patch] + gates[position:]}
                metrics = evaluate(trial)
                record = (metrics["off_diagonal_error"], position,
                          control, target, trial)
                if best is None or record[0] < best[0]:
                    best = record
                count += 1
    if args.output:
        args.output.write_text(json.dumps(best[4], indent=2) + "\n")
    print(json.dumps({"insertions_tested": count,
                      "best_off_diagonal_error": best[0],
                      "best_position": best[1],
                      "best_control": best[2],
                      "best_target": best[3]}, indent=2))


if __name__ == "__main__":
    main()
