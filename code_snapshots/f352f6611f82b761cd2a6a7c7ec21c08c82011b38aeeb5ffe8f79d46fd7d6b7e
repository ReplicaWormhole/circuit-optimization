"""Bounded two-entangler CX-power homotopy from exact14 toward 13 CX.

Drive exact14 CX6 from alpha=1 to alpha=0 via H CP(pi*alpha) H. At each
intermediate stage, independently optimize the power beta of CX7 or CX8
along with all local SU(2) layers, with a growing quadratic penalty that
returns beta toward one. At the final stage beta is fixed exactly at one,
so the serialized endpoint has thirteen actual CNOTs. This is a numerical
continuation experiment and does not establish a lower bound if it fails.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import embedded_one_qubit, evaluate, right_shift
from delete14_search import axis_rotation, fixed_matrix, kron_gate
from search13_cx_power_homotopy import (
    BASE, powered_cx, serialize_endpoint, split_exact14,
)


ROOT = Path(__file__).resolve().parent
REMOVED = 6
PARTNERS = (7, 8)
STAGES = ((0.75, 0.01), (0.5, 0.03), (0.35, 0.1),
          (0.2, 0.3), (0.0, None))
SIGMAS = (0.02, 0.15)
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64


def partner_data(gate):
    c, t = gate["control"], gate["target"]
    h = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
    ht = torch.as_tensor(embedded_one_qubit(h, t, 4), dtype=CDTYPE)
    p11 = torch.as_tensor([int(bool(basis & (8 >> c)) and
                               bool(basis & (8 >> t)))
                            for basis in range(16)], dtype=CDTYPE)
    return ht, p11


def partner_matrix(beta, ht, p11):
    phase = torch.exp((1j * math.pi) * beta)
    return ht @ torch.diag(1 + (phase - 1) * p11) @ ht


def fit_stage(angles0, beta0, alpha, penalty, partner, cxs, fixed, maxiter):
    free_beta = penalty is not None
    cxm = [torch.as_tensor(fixed_matrix(g), dtype=CDTYPE) for g in cxs]
    cxm[REMOVED] = powered_cx(cxs[REMOVED], alpha)
    ht, p11 = partner_data(cxs[partner])
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)

    def objective(flat):
        params = torch.tensor(flat, dtype=RDTYPE, requires_grad=True)
        angles = params[:180].reshape(15, 4, 3)
        beta = params[180] if free_beta else torch.as_tensor(1.0, dtype=RDTYPE)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(15):
            u = fixed[slot] @ u
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
            if slot < 14:
                gate = partner_matrix(beta, ht, p11) if slot == partner else cxm[slot]
                u = gate @ u
        d = u @ v @ u.conj().T
        off = d - torch.diag(torch.diagonal(d))
        raw_loss = (off.abs() ** 2).sum().real / 16
        total_loss = raw_loss + penalty * (beta - 1).square() if free_beta else raw_loss
        grad = torch.autograd.grad(total_loss, params)[0]
        return float(total_loss.detach()), grad.detach().numpy().copy()

    x0 = angles0.ravel()
    if free_beta:
        x0 = np.append(x0, beta0)
        bounds = [(None, None)] * 180 + [(0.0, 2.0)]
    else:
        bounds = None
    initial = objective(x0)[0]
    opt = minimize(objective, x0, method="L-BFGS-B", jac=True, bounds=bounds,
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    angles = opt.x[:180].reshape(15, 4, 3)
    beta = float(opt.x[180]) if free_beta else 1.0
    raw_loss = float(opt.fun - penalty * (beta - 1) ** 2) if free_beta else float(opt.fun)
    return angles, beta, {"alpha": alpha, "beta": beta,
                          "beta_penalty": penalty, "initial_total_loss": initial,
                          "raw_cycle_loss": raw_loss, "total_loss": float(opt.fun),
                          "nit": int(opt.nit), "nfev": int(opt.nfev),
                          "optimizer_success": bool(opt.success),
                          "optimizer_message": str(opt.message)}


def save_snapshot(prefix, partner, start, alpha, angles, beta):
    path = ROOT / f"{prefix}_p{partner}_s{start}_a{alpha:.3f}.npz"
    np.savez(path, angles=np.asarray(angles, dtype=np.float64),
             beta=np.asarray(beta, dtype=np.float64))
    return {"path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=14500)
    parser.add_argument("--maxiter", type=int, default=150)
    parser.add_argument("--prefix", default="search13_two_cx_power_homotopy")
    args = parser.parse_args()
    chunks, cxs, fixed = split_exact14()
    rows = []
    for partner_index, partner in enumerate(PARTNERS):
        for start, sigma in enumerate(SIGMAS):
            seed = args.seed + 2 * partner_index + start
            rng = np.random.default_rng(seed)
            angles = rng.normal(0, sigma, (15, 4, 3))
            beta = 1.0
            stages = []
            for alpha, penalty in STAGES:
                angles, beta, record = fit_stage(angles, beta, alpha, penalty,
                                                 partner, cxs, fixed, args.maxiter)
                record["snapshot"] = save_snapshot(args.prefix, partner, start,
                                                     alpha, angles, beta)
                stages.append(record)
                print(json.dumps({"partner": partner, "seed": seed,
                                  "stage": record}, sort_keys=True), flush=True)
            candidate = serialize_endpoint(chunks, cxs, angles, REMOVED)
            path = ROOT / f"{args.prefix}_p{partner}_s{start}_endpoint.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            check = evaluate(candidate)
            row = {"partner": partner, "seed": seed, "sigma": sigma,
                   "stages": stages, "candidate": path.name,
                   "endpoint_off_diagonal_error": check["off_diagonal_error"],
                   "endpoint_eigenvalue_root_error": check["eigenvalue_root_error"],
                   "endpoint_valid_diagonalizer": check["valid_diagonalizer"]}
            rows.append(row)
            print(json.dumps({"endpoint": row}, sort_keys=True), flush=True)
    result = {"base": BASE.name, "removed_cx": REMOVED,
              "partners": PARTNERS, "stages": STAGES, "sigmas": SIGMAS,
              "maxiter_per_stage": args.maxiter, "rows": rows,
              "best": min(rows, key=lambda r: r["stages"][-1]["raw_cycle_loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name,
                      "best_candidate": result["best"]["candidate"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
