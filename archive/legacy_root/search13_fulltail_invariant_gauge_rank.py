"""Bounded continuous output-eigenspace gauge rank search for exact14's tail.

This is a numerical screen, not an exact rank certificate. It optimizes the
sum of squares of the eight smallest balanced-cut realignment singular values
over the full block-unitary output gauge U(6)xU(4)xU(3)xU(3).
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import right_shift
from delete14_search import fixed_matrix

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_fulltail_invariant_gauge_rank_result.json"
CUTS = ((0, 1), (0, 2), (0, 3))
torch.set_num_threads(1)


def circuit_matrix(gates):
    out = np.eye(16, dtype=complex)
    for gate in gates:
        out = fixed_matrix(gate) @ out
    return out


def setup():
    gates = json.loads(SOURCE.read_text())["gates"]
    prefix = circuit_matrix(gates[:13])
    tail = circuit_matrix(gates[13:])
    full = tail @ prefix
    diagonal = np.diag(full @ right_shift(4) @ full.conj().T)
    roots = (1, 1j, -1, -1j)
    labels = [int(np.argmin([abs(z - root) for root in roots]))
              for z in diagonal]
    assert max(abs(z - roots[k]) for z, k in zip(diagonal, labels)) < 1e-12
    sectors = [tuple(j for j, label in enumerate(labels) if label == k)
               for k in range(4)]
    assert sorted(map(len, sectors)) == [3, 3, 4, 6]
    return tail, prefix, labels, sectors


def gauge(angles, sectors):
    out = torch.zeros((16, 16), dtype=torch.complex128)
    pos = 0
    for sector in sectors:
        n = len(sector)
        block = torch.zeros((n, n), dtype=torch.complex128)
        for i in range(n):
            block[i, i] = angles[pos]
            pos += 1
        for i in range(n):
            for j in range(i + 1, n):
                z = angles[pos] + 1j * angles[pos + 1]
                block[i, j] = z
                block[j, i] = torch.conj(z)
                pos += 2
        unitary = torch.matrix_exp(1j * block)
        indices = torch.as_tensor(sector, dtype=torch.long)
        out[indices[:, None], indices[None, :]] = unitary
    assert pos == angles.numel()
    return out


def realign(matrix, cut):
    left = tuple(cut) + tuple(4 + j for j in cut)
    right = tuple(j for j in range(8) if j not in left)
    return matrix.reshape((2,) * 8).permute(left + right).reshape(16, 16)


def run(cut, seed, sigma, maxiter, tail, prefix, labels, sectors):
    rng = np.random.default_rng(seed)
    nvar = sum(len(group) ** 2 for group in sectors)
    initial = rng.normal(0, sigma, nvar)
    target = torch.as_tensor(tail, dtype=torch.complex128)

    def objective(x):
        angles = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        gauged_tail = gauge(angles, sectors) @ target
        singular = torch.linalg.svdvals(realign(gauged_tail, cut))
        loss = (singular[8:] ** 2).sum().real / 16
        gradient = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), gradient.detach().numpy().copy()

    initial_loss = objective(initial)[0]
    fit = minimize(objective, initial, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-14,
                            "gtol": 1e-10, "maxls": 30})
    with torch.no_grad():
        angles = torch.as_tensor(fit.x, dtype=torch.float64)
        g = gauge(angles, sectors).numpy()
    trial = g @ tail
    singular = np.linalg.svd(realign(torch.as_tensor(trial), cut).numpy(),
                             compute_uv=False)
    whole = trial @ prefix
    off = whole @ right_shift(4) @ whole.conj().T
    off -= np.diag(np.diag(off))
    return {
        "cut": list(cut), "seed": seed, "sigma": sigma,
        "maxiter": maxiter, "initial_loss": initial_loss,
        "loss": float(fit.fun), "iterations": int(fit.nit),
        "optimizer_success": bool(fit.success), "optimizer_message": str(fit.message),
        "singular_values": singular.tolist(),
        "rank_tolerance_1e_8": int(np.count_nonzero(singular > 1e-8)),
        "smallest_singular": float(singular[-1]),
        "eighth_ninth_gap": float(singular[7] - singular[8]),
        "gauge_unitarity_error": float(np.max(np.abs(g @ g.conj().T - np.eye(16)))),
        "full_circuit_offdiagonal_error": float(np.max(np.abs(off))),
        "angles": fit.x.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--seed", type=int, default=25100)
    parser.add_argument("--sigma", type=float, default=0.3)
    args = parser.parse_args()
    tail, prefix, labels, sectors = setup()
    baseline = [np.linalg.svd(realign(torch.as_tensor(tail), cut).numpy(),
                              compute_uv=False).tolist() for cut in CUTS]
    rows = []
    for cut_index, cut in enumerate(CUTS):
        for start in range(args.starts):
            row = run(cut, args.seed + args.starts * cut_index + start,
                      args.sigma, args.maxiter, tail, prefix, labels, sectors)
            rows.append(row)
            print(json.dumps({key: row[key] for key in ("cut", "seed", "loss",
                  "rank_tolerance_1e_8", "smallest_singular", "optimizer_success")}),
                  flush=True)
    result = {
        "source": SOURCE.name, "prefix_cnot_count": 6,
        "tail_cnot_count": 8, "output_eigenvalue_labels": labels,
        "gauge_sector_dimensions": list(map(len, sectors)),
        "objective": "sum of squared eight smallest balanced realignment singular values / 16",
        "baseline_singular_values": baseline,
        "starts_per_cut": args.starts, "maxiter": args.maxiter,
        "rows": rows,
        "best_by_cut": [min((r for r in rows if r["cut"] == list(cut)),
                            key=lambda r: r["loss"]) for cut in CUTS],
        "scope": "bounded numerical rank screen; does not certify exact rank or synthesis",
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
