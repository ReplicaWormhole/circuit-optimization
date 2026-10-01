"""Bounded simultaneous balanced-cut rank search around the exact gauge.

Apply arbitrary block-unitary output gauges preserving the four shift
eigenvalue sectors to the exact rank-eight-gauged eight-CNOT tail.  Minimize
the squared singular values below rank eight across cuts 01|23 and 03|12.
This is a numerical screen; low loss alone would not certify exact rank.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import right_shift
from search13_fulltail_invariant_gauge_rank import gauge, realign, setup
from search13_rank8_bridge import exact_gauge


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_doublecut_gauge_fit_result.json"
CUTS = ((0, 1), (0, 3))
torch.set_num_threads(1)


def optimize(initial, target, sectors, weights, maxiter):
    def objective(x):
        angles = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        trial = gauge(angles, sectors) @ target
        losses = []
        for cut in CUTS:
            singular = torch.linalg.svdvals(realign(trial, cut))
            losses.append((singular[8:]**2).sum().real/16)
        loss = sum(w*l for w, l in zip(weights, losses))
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().copy()

    initial_loss = objective(initial)[0]
    fit = minimize(objective, initial, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    return fit.x, {"weights": list(weights), "start_loss": initial_loss,
                   "loss": float(fit.fun), "iterations": int(fit.nit),
                   "success": bool(fit.success), "message": str(fit.message)}


def inspect(angles, target, prefix, sectors):
    with torch.no_grad():
        h = gauge(torch.as_tensor(angles, dtype=torch.float64), sectors).numpy()
    trial = h @ target
    singular = {}
    losses = {}
    ranks = {}
    for cut in CUTS:
        key = "".join(map(str, cut))
        matrix = realign(torch.as_tensor(trial), cut).numpy()
        sv = np.linalg.svd(matrix, compute_uv=False)
        singular[key] = sv.tolist()
        losses[key] = float(np.sum(sv[8:]**2)/16)
        ranks[key] = int(np.count_nonzero(sv > 1e-8))
    whole = trial @ prefix
    transformed = whole @ right_shift(4) @ whole.conj().T
    off = transformed - np.diag(np.diag(transformed))
    return {"singular_values": singular, "rank_at_1e_minus_8": ranks,
            "tail_losses": losses,
            "total_loss": sum(losses.values()),
            "gauge_unitarity_error": float(np.max(np.abs(h @ h.conj().T - np.eye(16)))),
            "full_diagonalization_offdiag": float(np.max(np.abs(off)))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=29101)
    parser.add_argument("--stage1-maxiter", type=int, default=120)
    parser.add_argument("--stage2-maxiter", type=int, default=180)
    args = parser.parse_args()
    tail, prefix, labels, sectors = setup()
    g = exact_gauge()
    assert all(abs(g[i, j]) < 1e-14 or labels[i] == labels[j]
               for i in range(16) for j in range(16))
    base = g @ tail
    base_inspection = inspect(np.zeros(sum(len(s)**2 for s in sectors)),
                              base, prefix, sectors)
    rng = np.random.default_rng(args.seed)
    nvar = sum(len(s)**2 for s in sectors)
    starts = (("exact", np.zeros(nvar)),
              ("near", rng.normal(0, 0.06, nvar)),
              ("wide", rng.normal(0, 0.2, nvar)))
    target = torch.as_tensor(base, dtype=torch.complex128)
    rows = []
    for name, initial in starts:
        angles, first = optimize(initial, target, sectors, (1.0, 4.0),
                                  args.stage1_maxiter)
        angles, second = optimize(angles, target, sectors, (1.0, 1.0),
                                   args.stage2_maxiter)
        observation = inspect(angles, base, prefix, sectors)
        row = {"start": name, "stages": [first, second],
               "angles": angles.tolist(), **observation}
        rows.append(row)
        print(json.dumps({k: row[k] for k in ("start", "tail_losses",
                        "rank_at_1e_minus_8", "total_loss",
                        "full_diagonalization_offdiag")}, sort_keys=True), flush=True)
    result = {"source": "topology14_exact_matchgate_rational.json",
              "exact_gauge_source": "search13_fulltail_invariant_exact_witness.py",
              "seed": args.seed, "stage1_maxiter": args.stage1_maxiter,
              "stage2_maxiter": args.stage2_maxiter,
              "sector_dimensions": list(map(len, sectors)),
              "cuts": [list(c) for c in CUTS],
              "objective": "sum of squared singular values 9..16 / 16 on each cut",
              "baseline": base_inspection, "rows": rows,
              "best_start": min(rows, key=lambda r: r["total_loss"])["start"],
              "scope": "three bounded numerical output-gauge fits; not an exact rank or CNOT certificate"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
