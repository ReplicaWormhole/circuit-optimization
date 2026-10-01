"""Rescaled numerical refinement of the simultaneous mixed-cut gauge lead.

The run293 joint-exact start has rank-eight tail losses near 2e-15 at
masks 83 and 163. Multiplying the objective and gradient by a large
constant prevents a small absolute objective from prematurely ending
L-BFGS-B. This is numerical evidence, not an exact rank certificate.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from search13_fulltail_invariant_gauge_rank import gauge, setup
from search13_mixedcut_gauge_rank import MASKS, inspect, loss_by_mask
from search13_rank8_bridge import exact_gauge


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_mixedcut_gauge_rank_result.json"
OUT = ROOT / "search13_mixedcut_gauge_refine_result.json"
torch.set_num_threads(1)


def optimize(initial, target, sectors, scale, maxiter):
    target = torch.as_tensor(target, dtype=torch.complex128)

    def objective(x):
        angles = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        losses = loss_by_mask(angles, target, sectors)
        value = scale * sum(losses.values())
        gradient = torch.autograd.grad(value, angles)[0]
        return float(value.detach()), gradient.detach().numpy().copy()

    fit = minimize(objective, initial, method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-16,
                            "gtol": 1e-10, "maxls": 40})
    return fit.x, {"scale": scale, "scaled_loss": float(fit.fun),
                   "iterations": int(fit.nit), "message": str(fit.message),
                   "success": bool(fit.success)}


def all_mixed_ranks(matrix, threshold):
    tensor = matrix.reshape((2,) * 8)
    rows = {}
    for mask in range(1, 255, 2):
        left = tuple(j for j in range(8) if mask & (1 << j))
        right = tuple(j for j in range(8) if j not in left)
        flat = tensor.transpose(left + right).reshape(1 << len(left), 1 << len(right))
        sv = np.linalg.svd(flat, compute_uv=False)
        rows[str(mask)] = int(np.count_nonzero(sv > threshold))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=250)
    args = parser.parse_args()
    previous = json.loads(SOURCE.read_text())
    start = next(row for row in previous["rows"]
                 if row["goal"] == "joint" and row["start"] == "exact")
    angles = np.array(start["angles"], dtype=float)
    tail, prefix, labels, sectors = setup()
    g = exact_gauge()
    assert all(abs(g[i, j]) < 1e-14 or labels[i] == labels[j]
               for i in range(16) for j in range(16))
    base = g @ tail
    stages = []
    for scale in (1e8, 1e12, 1e16):
        angles, fit = optimize(angles, base, sectors, scale, args.maxiter)
        stages.append({"fit": fit, "observation": inspect(angles, base, prefix, sectors)})
        print(json.dumps({"scale": scale, "scaled_loss": fit["scaled_loss"],
                          "cuts": {k: v["loss_below_rank8"]
                                   for k, v in stages[-1]["observation"]["cuts"].items()}},
                         sort_keys=True), flush=True)
    with torch.no_grad():
        refined = gauge(torch.as_tensor(angles, dtype=torch.float64), sectors).numpy() @ base
    result = {"source": SOURCE.name, "masks": MASKS,
              "maxiter_per_stage": args.maxiter, "stages": stages,
              "refined_angles": angles.tolist(),
              "numerical_rank_profiles_at_1e_minus_6": all_mixed_ranks(refined, 1e-6),
              "numerical_rank_profiles_at_1e_minus_9": all_mixed_ranks(refined, 1e-9),
              "scope": "numerical refinement only; rank profiles depend on thresholds"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
