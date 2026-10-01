"""Large-noise multistarts on the two lowest final-loss run335 topologies."""

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
SOURCE = ROOT / "search13_joint_prefix_slot8_refine_trial1.json"
RESULT = ROOT / "search13_joint_prefix_two_slot_v2_result.json"
FROZEN = ROOT / "search13_joint_prefix_two_slot_v2_frozen_proposals.json"
OUT = ROOT / "search13_joint_prefix_two_slot_multistart_result.json"
torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=33700)
    parser.add_argument("--maxiter", type=int, default=400)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    angles = decode_layers(source)
    prior = json.loads(RESULT.read_text())
    selected = sorted(prior["rows"], key=lambda row: row["loss"])[:2]
    assert [row["case"] for row in selected] == [15, 9]
    frozen = {row["case"]: row["schedule"]
              for row in json.loads(FROZEN.read_text())["proposals"]}
    deleted, rewire, blocks = CASES[1]
    chunks, cxs, fixed, _ = prepare(deleted, blocks, prefix_rewire=rewire)
    original_schedule = [(g["control"], g["target"]) for g in cxs]
    original_schedule[8] = (1, 2)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed,
                          matrix_for_schedule(original_schedule), v, False)
    assert abs(roundtrip - 0.06019279505817157) < 1e-10, roundtrip
    rows = []
    for rank, chosen in enumerate(selected):
        schedule = tuple(map(tuple, frozen[chosen["case"]]))
        for start in range(2):
            seed = args.seed + 100 * rank + start
            initial = angles + np.random.default_rng(seed).normal(0, 0.8, angles.shape)
            fitted, record = fit(schedule, initial, fixed, v, args.maxiter)
            candidate = serialize(chunks, schedule, fitted)
            path = ROOT / f"search13_joint_prefix_two_slot_multistart_case{chosen['case']}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            check = evaluate(candidate)
            row = {"case": chosen["case"], "start": start, "seed": seed,
                   "noise_sigma": 0.8, "schedule": schedule, **record,
                   "candidate": path.name, "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid_diagonalizer": check["valid_diagonalizer"]}
            rows.append(row)
            print(json.dumps({"case": row["case"], "start": start,
                              "initial_loss": row["initial_loss"],
                              "loss": row["loss"],
                              "off_diagonal_error": row["off_diagonal_error"],
                              "valid": row["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["off_diagonal_error"])
    result = {"source": SOURCE.name, "topology_selection": "two lowest final fitted Frobenius losses from run 335, cases 15 and 9",
              "frozen_proposals": FROZEN.name,
              "initialization": "run 332 source optimized local angles plus independent Gaussian noise of std 0.8; not small perturbation or fitted run 335 locals",
              "source_roundtrip_loss": roundtrip,
              "seed": args.seed, "maxiter": args.maxiter, "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "four numerical multistarts on two fixed topologies, full 168 local parameter freedom; failure is not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
