"""Deepen the best changed-prefix near-hit from run 325.

Continue the 13-CNOT case-1 random-start candidate without perturbation and
with two seeded perturbations. Every local gate in the circuit is refitted.
This is a bounded numerical refinement and cannot establish a no-go result.
"""

import argparse
import json
from pathlib import Path

from check_circuit import evaluate
from search13_joint_prefix_orbit_fit import CASES
from search13_multi_pair_block_rewire import fit
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_joint_prefix_orbit_case1_random_direct.json"
PRIOR = ROOT / "search13_joint_prefix_orbit_fit_result.json"
OUT = ROOT / "search13_joint_prefix_orbit_refine_result.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=32900)
    parser.add_argument("--maxiter", type=int, default=800)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    prior = json.loads(PRIOR.read_text())
    prior_row = next(r for r in prior["rows"]
                     if r["case"] == 1 and r["start"] == "random")
    expected = prior_row["direct"]["loss"]
    angles = decode_layers(source)
    deleted, rewire, blocks = CASES[1]
    rows = []
    for trial, sigma in enumerate((0.0, 0.02, 0.08)):
        candidate, record = fit(
            deleted, blocks, args.seed + trial, sigma, args.maxiter,
            initial_angles=angles, prefix_rewire=rewire)
        if trial == 0:
            assert abs(record["initial_loss"] - expected) < 1e-8, (
                record["initial_loss"], expected)
        path = ROOT / f"search13_joint_prefix_orbit_refine_{trial}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"trial": trial, "sigma": sigma,
               "candidate": path.name, **record,
               "cnot_count": check["cnot_count"]}
        rows.append(row)
        print(json.dumps({"trial": trial, "sigma": sigma,
                          "initial_loss": record["initial_loss"],
                          "loss": record["loss"],
                          "off_diagonal_error": check["off_diagonal_error"],
                          "valid": check["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda r: r["off_diagonal_error"])
    result = {"source": SOURCE.name, "seed": args.seed,
              "maxiter": args.maxiter, "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "one whole-circuit 13-CNOT topology, three numerical continuations; failure is not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
