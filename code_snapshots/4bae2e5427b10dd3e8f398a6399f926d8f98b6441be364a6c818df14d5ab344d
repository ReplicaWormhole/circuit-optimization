"""Insert one CZ into the best 15-CNOT Bell-replacement candidate."""

import argparse
import json
from pathlib import Path

from check_circuit import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).parent
    candidate = json.loads((root / "last_f2_swap_basis_best.json").read_text())
    assert sum(gate["gate"] == "cx" for gate in candidate["gates"]) == 15
    gates = candidate["gates"]
    best = None
    tested = 0
    for position in range(len(gates) + 1):
        for first in range(4):
            for second in range(first + 1, 4):
                patch = [{"gate": "h", "qubit": second},
                         {"gate": "cx", "control": first, "target": second},
                         {"gate": "h", "qubit": second}]
                trial = {"n": 4, "gates": gates[:position] + patch + gates[position:]}
                metrics = evaluate(trial)
                record = (metrics["off_diagonal_error"], position,
                          (first, second), trial)
                if best is None or record[0] < best[0]:
                    best = record
                tested += 1
    if args.output:
        args.output.write_text(json.dumps(best[3], indent=2) + "\n")
    print(json.dumps({"insertions_tested": tested,
                      "best_off_diagonal_error": best[0],
                      "best_position": best[1],
                      "best_pair": best[2]}, indent=2))


if __name__ == "__main__":
    main()
