"""Positive output-row residual weights followed by ordinary cycle polish.

Eight trials on two fixed 13-CNOT topologies. All 168 local SU(2)
parameters are free. Four fixed positive weight vectors per topology:
source-residual adaptive, reversed emphasis, and two seeded random vectors.
The weighted objective has the same exact zeros as the ordinary loss.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import matrix_for_schedule, objective, serialize
from search13_fulltail_invariant_gauge_rank import circuit_matrix
from search13_joint_prefix_orbit_fit import CASES, fit, unitary
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_positive_row_weights_result.json"
SOURCES = ("search13_joint_prefix_slot8_refine_trial1.json",
           "search13_joint_prefix_two_slot_v2_case15.json")
FROZEN = ROOT / "search13_joint_prefix_two_slot_v2_frozen_proposals.json"
torch.set_num_threads(1)


def weighted_objective(flat, fixed, cxm, v, weights):
    angles = torch.tensor(np.asarray(flat).reshape(14, 4, 3),
                          dtype=torch.float64, requires_grad=True)
    u = unitary(angles, fixed, cxm)
    d = u @ v @ u.conj().T
    off = d - torch.diag(torch.diagonal(d))
    loss = (weights * (off.abs() ** 2).sum(dim=1)).sum().real / weights.sum()
    grad = torch.autograd.grad(loss, angles)[0]
    return float(loss.detach()), grad.detach().numpy().ravel().copy()


def weight_vectors(candidate, seed):
    u = circuit_matrix(candidate["gates"])
    d = u @ right_shift(4) @ u.conj().T
    off = d - np.diag(np.diag(d))
    residual = np.sum(np.abs(off) ** 2, axis=1)
    assert max(residual) > 1e-16
    adaptive = 1 + 9 * residual / max(residual)
    rows = [("residual_adaptive", None, adaptive),
            ("reversed_emphasis", None, 11 - adaptive)]
    for start in range(2):
        rng = np.random.default_rng(seed + start)
        raw = np.exp(rng.normal(0, 1, 16))
        weights = np.minimum(raw / min(raw), 10)
        rows.append(("positive_random", seed + start, weights))
    return residual, rows


def checked(path, candidate):
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    q = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": q["cnot_count"],
            "off_diagonal_error": q["off_diagonal_error"],
            "valid_diagonalizer": q["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=33900)
    parser.add_argument("--weighted-maxiter", type=int, default=200)
    parser.add_argument("--direct-maxiter", type=int, default=250)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    deleted, rewire, blocks = CASES[1]
    chunks, old_cxs, fixed, _ = prepare(deleted, blocks, prefix_rewire=rewire)
    base_schedule = [(g["control"], g["target"]) for g in old_cxs]
    base_schedule[8] = (1, 2)
    case15 = next(row["schedule"] for row in json.loads(FROZEN.read_text())["proposals"]
                  if row["case"] == 15)
    schedules = (base_schedule, tuple(map(tuple, case15)))
    expected_losses = (0.06019279505817157, 0.060192795089502246)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    rows = []
    sources = []
    for topology, (source_name, schedule, expected) in enumerate(
            zip(SOURCES, schedules, expected_losses)):
        path = ROOT / source_name
        candidate = json.loads(path.read_text())
        assert tuple(schedule) == tuple((g["control"], g["target"])
                                        for g in candidate["gates"] if g["gate"] == "cx")
        angles = decode_layers(candidate)
        cxm = matrix_for_schedule(schedule)
        direct_loss, direct_grad = objective(angles.ravel(), fixed, cxm, v, True)
        assert abs(direct_loss - expected) < 1e-10, (direct_loss, expected)
        ones_loss, ones_grad = weighted_objective(
            angles.ravel(), fixed, cxm, v, torch.ones(16, dtype=torch.float64))
        assert abs(ones_loss - direct_loss) < 1e-12
        assert np.max(np.abs(ones_grad - direct_grad)) < 1e-12
        residual, vectors = weight_vectors(candidate, args.seed + 100 * topology)
        sources.append({"source": source_name,
                        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "schedule": schedule,
                        "roundtrip_loss": direct_loss,
                        "roundtrip_error": abs(direct_loss - expected),
                        "source_row_residuals": residual.tolist()})
        print(json.dumps({"topology": topology, "source": source_name,
                          "roundtrip_loss": direct_loss,
                          "roundtrip_error": abs(direct_loss - expected)}), flush=True)
        for trial, (kind, seed, weights_np) in enumerate(vectors):
            assert min(weights_np) > 0
            assert max(weights_np) / min(weights_np) <= 10 + 1e-12
            weights = torch.as_tensor(weights_np, dtype=torch.float64)
            weighted_angles, weighted = fit(
                lambda x: weighted_objective(x, fixed, cxm, v, weights),
                angles, args.weighted_maxiter)
            weighted_candidate = serialize(chunks, schedule, weighted_angles)
            weighted_path = ROOT / f"search13_positive_row_weights_t{topology}_w{trial}_weighted.json"
            weighted_check = checked(weighted_path, weighted_candidate)
            direct_angles, direct = fit(
                lambda x: objective(x, fixed, cxm, v, True),
                weighted_angles, args.direct_maxiter)
            direct_candidate = serialize(chunks, schedule, direct_angles)
            direct_path = ROOT / f"search13_positive_row_weights_t{topology}_w{trial}_direct.json"
            direct_check = checked(direct_path, direct_candidate)
            row = {"topology": topology, "trial": trial, "weight_kind": kind,
                   "random_weight_seed": seed, "weights": weights_np.tolist(),
                   "weighted": {**weighted, **weighted_check},
                   "direct": {**direct, **direct_check}}
            rows.append(row)
            print(json.dumps({"topology": topology, "trial": trial,
                              "weight_kind": kind,
                              "weighted_loss": weighted["loss"],
                              "direct_loss": direct["loss"],
                              "off_diagonal_error": direct_check["off_diagonal_error"],
                              "valid": direct_check["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["direct"]["off_diagonal_error"])
    result = {"seed": args.seed, "sources": sources,
              "weight_rule": "1+9*row residual/max residual, reverse11-w, and clipped exp-normal vectors with minimum1 and maximum10",
              "objective": "sum_i w_i sum_j!=i |(U V4 Udagger)_ij|^2 / sum_i w_i",
              "strictly_positive_weights": True, "fixed_eigenvalue_labels": False,
              "weighted_maxiter": args.weighted_maxiter,
              "direct_maxiter": args.direct_maxiter, "rows": rows,
              "best_candidate": best["direct"]["candidate"],
              "best_off_diagonal_error": best["direct"]["off_diagonal_error"],
              "scope": "eight fixed positive-row-weight trials on two 13-CNOT topologies; full local freedom; numerical failures are not exclusions"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
