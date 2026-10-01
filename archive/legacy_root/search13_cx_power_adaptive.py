"""Finer CX-power continuation from run214's near-exact alpha=0.5 branch.

Deterministically replay position 6, seed 14202 through alpha=0.5 using
run214's 100-iteration stages, then preserve that exact NumPy parameter
snapshot. Continue to zero in smaller alpha steps with 250-iteration
SU(2)-layer fits. At alpha=0.35 and 0.25, test two deterministic small
perturbations as well as the unperturbed warm start and keep the best.
This is a bounded numerical path, not a proof of existence/nonexistence.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from search13_cx_power_homotopy import (
    BASE, fit_stage, serialize_endpoint, split_exact14,
)


ROOT = Path(__file__).resolve().parent
POSITION = 6
REPLAY_SEED = 14202
REPLAY_SIGMA = 0.02
REPLAY_STAGES = (0.875, 0.75, 0.5)
CONTINUATION_STAGES = (0.45, 0.4, 0.35, 0.3, 0.25,
                       0.2, 0.15, 0.1, 0.05, 0.0)
PERTURB_STAGES = (0.35, 0.25)
PERTURB_SIGMAS = (0.01, 0.04)


def save_angles(prefix, alpha, angles):
    path = ROOT / f"{prefix}_alpha{alpha:.3f}_angles.npy"
    np.save(path, np.asarray(angles, dtype=np.float64), allow_pickle=False)
    raw = path.read_bytes()
    return {"path": path.name, "sha256": hashlib.sha256(raw).hexdigest(),
            "shape": list(angles.shape)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=250)
    parser.add_argument("--perturb-seed", type=int, default=14300)
    parser.add_argument("--prefix", default="search13_cx_power_adaptive")
    args = parser.parse_args()
    chunks, cxs, fixed = split_exact14()
    rng = np.random.default_rng(REPLAY_SEED)
    angles = rng.normal(0, REPLAY_SIGMA, (15, 4, 3))
    replay = []
    for alpha in REPLAY_STAGES:
        angles, row = fit_stage(angles, alpha, POSITION, cxs, fixed, 100)
        replay.append(row)
        print(json.dumps({"replay": row}, sort_keys=True), flush=True)
    previous = json.loads((ROOT / "search13_cx_power_homotopy_run214_result.json").read_text())
    expected = next(r for r in previous["rows"] if r["removed"] == POSITION
                    and r["seed"] == REPLAY_SEED)["stages"][:3]
    replay_error = max(abs(r["loss"] - e["loss"])
                       for r, e in zip(replay, expected))
    if replay_error > 1e-8:
        raise AssertionError(f"deterministic replay drifted by {replay_error}")
    snapshots = {"0.500": save_angles(args.prefix, 0.5, angles)}
    continuation = []
    perturb_rng = np.random.default_rng(args.perturb_seed)
    for alpha in CONTINUATION_STAGES:
        starts = [("warm", angles)]
        if alpha in PERTURB_STAGES:
            starts.extend((f"perturb_{sigma}",
                           angles + perturb_rng.normal(0, sigma, angles.shape))
                          for sigma in PERTURB_SIGMAS)
        trials = []
        for label, x0 in starts:
            candidate_angles, row = fit_stage(x0, alpha, POSITION, cxs,
                                               fixed, args.maxiter)
            trials.append((candidate_angles, {"start": label, **row}))
            print(json.dumps({"alpha": alpha, "trial": trials[-1][1]},
                             sort_keys=True), flush=True)
        angles, best = min(trials, key=lambda pair: pair[1]["loss"])
        snapshot = save_angles(args.prefix, alpha, angles)
        snapshots[f"{alpha:.3f}"] = snapshot
        continuation.append({"alpha": alpha, "trials": [r for _, r in trials],
                             "selected_start": best["start"],
                             "selected_loss": best["loss"],
                             "snapshot": snapshot})
    endpoint = serialize_endpoint(chunks, cxs, angles, POSITION)
    candidate_path = ROOT / f"{args.prefix}_endpoint.json"
    candidate_path.write_text(json.dumps(endpoint, indent=2) + "\n")
    check = evaluate(endpoint)
    result = {"base": BASE.name, "removed_cx": POSITION,
              "replay_seed": REPLAY_SEED, "replay_sigma": REPLAY_SIGMA,
              "replay_maxiter": 100, "replay": replay,
              "replay_max_loss_difference": replay_error,
              "continuation_maxiter": args.maxiter,
              "perturb_seed": args.perturb_seed,
              "perturb_sigmas": PERTURB_SIGMAS,
              "continuation": continuation, "snapshots": snapshots,
              "endpoint_candidate": candidate_path.name,
              "endpoint_off_diagonal_error": check["off_diagonal_error"],
              "endpoint_eigenvalue_root_error": check[
                  "eigenvalue_root_error"],
              "endpoint_valid_diagonalizer": check["valid_diagonalizer"]}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name,
                      "endpoint_candidate": candidate_path.name,
                      "endpoint_loss": continuation[-1]["selected_loss"],
                      "endpoint_off_diagonal_error": check[
                          "off_diagonal_error"],
                      "endpoint_valid_diagonalizer": check[
                          "valid_diagonalizer"]}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
