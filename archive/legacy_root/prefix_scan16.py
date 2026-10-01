"""Check whether a prefix of the 17-CNOT circuit already diagonalizes V4."""

import json
import argparse
from pathlib import Path

from check_circuit import evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).parent
    circuit = json.loads((root / "algebraic_boundary_candidate.json").read_text())
    results = []
    for length in range(len(circuit["gates"])):
        candidate = {"n": 4, "gates": circuit["gates"][:length]}
        result = evaluate(candidate)
        if result["cnot_count"] < 17:
            results.append({"gate_count": length,
                            "cnot_count": result["cnot_count"],
                            "off_diagonal_error": result["off_diagonal_error"]})
    best = min(results, key=lambda result: result["off_diagonal_error"])
    if args.output:
        best_candidate = {"n": 4, "gates": circuit["gates"][:best["gate_count"]]}
        args.output.write_text(json.dumps(best_candidate, indent=2) + "\n")
    print(json.dumps({"prefixes_checked": len(results), "best": best}, indent=2))


if __name__ == "__main__":
    main()
