"""Deep whole-circuit refinement of run330's slot8 topology lead."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_joint_prefix_orbit_fit import CASES
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_joint_prefix_slot_diversity_case8.json"
PRIOR = ROOT / "search13_joint_prefix_slot_diversity_result.json"
OUT = ROOT / "search13_joint_prefix_slot8_refine_result.json"
torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=33400)
    parser.add_argument("--maxiter", type=int, default=1000)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    prior = json.loads(PRIOR.read_text())
    expected = next(row["loss"] for row in prior["rows"] if row["case"] == 8)
    angles = decode_layers(source)
    deleted, rewire, blocks = CASES[1]
    chunks, cxs, fixed, _ = prepare(deleted, blocks, prefix_rewire=rewire)
    schedule = [(g["control"], g["target"]) for g in cxs]
    schedule[8] = (1, 2)
    assert tuple(schedule) == tuple((g["control"], g["target"])
                                    for g in source["gates"] if g["gate"] == "cx")
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed,
                          matrix_for_schedule(schedule), v, False)
    assert abs(roundtrip - expected) < 1e-10, (roundtrip, expected)
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "source_expected_loss": expected,
                      "roundtrip_error": abs(roundtrip - expected)}), flush=True)
    rows = []
    for trial, sigma in enumerate((0.0, 0.02, 0.08, 0.2)):
        rng = np.random.default_rng(args.seed + trial)
        initial = angles + rng.normal(0, sigma, angles.shape)
        fitted, record = fit(schedule, initial, fixed, v, args.maxiter)
        candidate = serialize(chunks, schedule, fitted)
        path = ROOT / f"search13_joint_prefix_slot8_refine_trial{trial}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"trial": trial, "sigma": sigma, "seed": args.seed + trial,
               **record, "candidate": path.name,
               "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps({"trial": trial, "sigma": sigma,
                          "initial_loss": row["initial_loss"],
                          "loss": row["loss"], "nit": row["nit"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda r: r["off_diagonal_error"])
    result = {"source": SOURCE.name, "seed": args.seed,
              "maxiter": args.maxiter, "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "schedule": schedule, "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "one13-CNOT topology, four whole-circuit numerical continuations; not an exact circuit or no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
