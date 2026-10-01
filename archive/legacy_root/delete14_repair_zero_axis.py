"""Repair legacy tail-search candidates whose zero SU2 slots became Z.

Only the twenty frozen prefix correction gates are changed.  The original
candidate remains archived in the ledger for provenance.
"""

import argparse
import json
import math
from pathlib import Path

from check_circuit import evaluate


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    candidate = json.loads(args.source.read_text())
    count = 0
    for gate in candidate["gates"]:
        if (gate["gate"] == "u3" and gate["theta"] == 0 and
                gate["phi"] == 0 and abs(gate["lam"] + math.pi) < 1e-12):
            gate["lam"] = 0.0
            count += 1
    if count != 20:
        raise ValueError(f"expected twenty frozen zero-axis slots, found {count}")
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    result = evaluate(candidate)
    print(json.dumps({"source": str(args.source), "output": str(args.output),
                      "repaired_slots": count,
                      "off_diagonal_error": result["off_diagonal_error"],
                      "valid": result["valid_diagonalizer"]}, indent=2))


if __name__ == "__main__":
    main()
