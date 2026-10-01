"""Perturb and refit a near-hit from the prefix deletion screen.

The saved u3 layers are decoded back to SU(2) axis-angle coordinates.
Before perturbation, the restored loss must match the candidate's checked
loss; this guards against serialization or decomposition mistakes.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from delete14_search import u3_to_axis
from search13_multi_pair_block_rewire import fit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_prefix_delete_fit_p4_s1.json"
BLOCKS = ((0, 2), (1, 3), (0, 1), (2, 3))


def decode_layers(candidate):
    chunks = [[]]
    for gate in candidate["gates"]:
        if gate["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if len(chunks) != 14:
        raise AssertionError(f"expected fourteen local chunks, got {len(chunks)}")
    layers = []
    for chunk in chunks:
        rotations = chunk[-4:]
        if [g["gate"] for g in rotations] != ["u3"] * 4:
            raise AssertionError("candidate does not contain the expected u3 layer")
        if [g["qubit"] for g in rotations] != list(range(4)):
            raise AssertionError("unexpected qubit order")
        layers.append([u3_to_axis(gate) for gate in rotations])
    return np.asarray(layers, dtype=float)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--seed", type=int, default=134300)
    parser.add_argument("--prefix", default="search13_prefix_perturb")
    args = parser.parse_args()
    source = json.loads(SOURCE.read_text())
    layers = decode_layers(source)
    checks = evaluate(source)
    rows = []
    for start, sigma in enumerate((0.0, 0.01, 0.04, 0.1)):
        seed = args.seed + start
        candidate, record = fit(4, BLOCKS, seed, sigma, args.maxiter,
                                initial_angles=layers)
        if start == 0:
            # The fit objective is squared Frobenius off-diagonal norm /16;
            # the run-230 saved record is our independent initial target.
            reference = json.loads((ROOT / "search13_prefix_delete_fit_result.json")
                                   .read_text())
            expected = next(r["loss"] for r in reference["rows"]
                            if r["deleted"] == 4 and r["start"] == 1)
            if abs(record["initial_loss"] - expected) > 1e-8:
                raise AssertionError(("roundtrip loss mismatch",
                                      record["initial_loss"], expected))
        path = ROOT / f"{args.prefix}_s{start}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        record.update(start=start, candidate=path.name,
                      source_off_diagonal_error=checks["off_diagonal_error"])
        rows.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
    result = {"source": SOURCE.name, "seed": args.seed,
              "maxiter": args.maxiter, "perturb_sigmas": [0, 0.01, 0.04, 0.1],
              "rows": rows, "best": min(rows, key=lambda row: row["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": path.name,
                      "best": result["best"]["candidate"]}), flush=True)


if __name__ == "__main__":
    main()
