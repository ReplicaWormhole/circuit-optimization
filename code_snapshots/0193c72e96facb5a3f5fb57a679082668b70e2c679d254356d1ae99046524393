"""Bounded output-gauge rank search at two mixed input/output cuts.

The exact run266 gauge gives rank 14 on masks 83 and 163. We vary a
block-unitary within every shift eigenspace and minimize squared singular
values below rank eight at each mask, then jointly. Results are numerical
leads only; they cannot prove a gauge-independent lower bound.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import right_shift
from search13_fulltail_invariant_gauge_rank import gauge, setup
from search13_rank8_bridge import exact_gauge


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_mixedcut_gauge_rank_result.json"
MASKS = (83, 163)
torch.set_num_threads(1)


def realign(matrix, mask):
    left = tuple(j for j in range(8) if mask & (1 << j))
    right = tuple(j for j in range(8) if j not in left)
    return matrix.reshape((2,) * 8).permute(left + right).reshape(16, 16)


def loss_by_mask(angles, base, sectors):
    trial = gauge(angles, sectors) @ base
    result = {}
    for mask in MASKS:
        singular = torch.linalg.svdvals(realign(trial, mask))
        result[mask] = (singular[8:]**2).sum().real / 16
    return result


def optimize(initial, base, sectors, weights, maxiter):
    target = torch.as_tensor(base, dtype=torch.complex128)

    def objective(x):
        angles = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        losses = loss_by_mask(angles, target, sectors)
        value = sum(weights[mask] * losses[mask] for mask in MASKS)
        grad = torch.autograd.grad(value, angles)[0]
        return float(value.detach()), grad.detach().numpy().copy()

    start = objective(initial)[0]
    fit = minimize(objective, initial, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    return fit.x, {"weights": weights, "start_loss": start,
                   "final_loss": float(fit.fun), "iterations": int(fit.nit),
                   "success": bool(fit.success), "message": str(fit.message)}


def inspect(angles, base, prefix, sectors):
    with torch.no_grad():
        h = gauge(torch.as_tensor(angles, dtype=torch.float64), sectors).numpy()
    trial = h @ base
    rows = {}
    for mask in MASKS:
        axes = tuple(j for j in range(8) if mask & (1 << j))
        other = tuple(j for j in range(8) if j not in axes)
        flat = trial.reshape((2,) * 8).transpose(axes + other).reshape(16, 16)
        sv = np.linalg.svd(flat, compute_uv=False)
        rows[str(mask)] = {"loss_below_rank8": float(np.sum(sv[8:]**2)/16),
                           "rank_at_1e_minus_8": int(np.count_nonzero(sv > 1e-8)),
                           "ninth_singular": float(sv[8]),
                           "singular_values": sv.tolist()}
    whole = trial @ prefix
    d = whole @ right_shift(4) @ whole.conj().T
    off = d - np.diag(np.diag(d))
    return {"cuts": rows,
            "max_diagonalization_offdiag": float(np.max(np.abs(off))),
            "gauge_unitarity_error": float(np.max(np.abs(h @ h.conj().T - np.eye(16))))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=29301)
    parser.add_argument("--maxiter", type=int, default=160)
    args = parser.parse_args()
    tail, prefix, labels, sectors = setup()
    g = exact_gauge()
    assert all(abs(g[i, j]) < 1e-14 or labels[i] == labels[j]
               for i in range(16) for j in range(16))
    base = g @ tail
    nvar = sum(len(s)**2 for s in sectors)
    rng = np.random.default_rng(args.seed)
    starts = (("exact", np.zeros(nvar)),
              ("perturbed", rng.normal(0, 0.12, nvar)))
    goals = (("mask83", {83: 1.0, 163: 0.0}),
             ("mask163", {83: 0.0, 163: 1.0}),
             ("joint", {83: 1.0, 163: 1.0}))
    rows = []
    for goal, weights in goals:
        for name, initial in starts:
            angles, fit = optimize(initial, base, sectors, weights, args.maxiter)
            observation = inspect(angles, base, prefix, sectors)
            row = {"goal": goal, "start": name, "fit": fit,
                   "angles": angles.tolist(), **observation}
            rows.append(row)
            print(json.dumps({"goal": goal, "start": name,
                              "losses": {k: v["loss_below_rank8"]
                                         for k, v in observation["cuts"].items()},
                              "offdiag": observation["max_diagonalization_offdiag"]},
                             sort_keys=True), flush=True)
    result = {"source": "exact run266 output gauge times exact14 tail",
              "masks": MASKS, "sector_dimensions": list(map(len, sectors)),
              "seed": args.seed, "maxiter": args.maxiter,
              "baseline": inspect(np.zeros(nvar), base, prefix, sectors),
              "rows": rows,
              "scope": "bounded numerical search for output gauge rank reduction; no exact rank or CNOT certificate"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
