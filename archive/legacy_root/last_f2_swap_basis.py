"""Test replacing final fermionic two-mode Fourier gates by one-CNOT bases."""

import argparse
import itertools
import json
from pathlib import Path

from check_circuit import evaluate


def alternatives(first, second, original):
    yield "original", original
    for control, target in ((first, second), (second, first)):
        for hadamard in (first, second):
            yield f"cx{control}{target}-h{hadamard}", [
                {"gate": "cx", "control": control, "target": target},
                {"gate": "h", "qubit": hadamard},
            ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidate = json.loads((Path(__file__).parent / "algebraic_boundary_candidate.json").read_text())
    gates = candidate["gates"]
    head, first, second = gates[:39], gates[39:49], gates[49:59]
    assert len(gates) == 59
    assert sum(g["gate"] == "cx" for g in first) == 2
    assert sum(g["gate"] == "cx" for g in second) == 2
    results = []
    for (name_a, block_a), (name_b, block_b) in itertools.product(
            alternatives(0, 1, first), alternatives(2, 3, second)):
        if name_a == name_b == "original":
            continue
        trial = {"n": 4, "gates": head + block_a + block_b}
        metrics = evaluate(trial)
        results.append((metrics["off_diagonal_error"], name_a, name_b,
                        metrics["cnot_count"], trial))
    best = min(results, key=lambda item: item[0])
    if args.output:
        args.output.write_text(json.dumps(best[4], indent=2) + "\n")
    print(json.dumps({"variants_tested": len(results), "best": {
        "off_diagonal_error": best[0], "first": best[1],
        "second": best[2], "cnot_count": best[3]}}, indent=2))


if __name__ == "__main__":
    main()
