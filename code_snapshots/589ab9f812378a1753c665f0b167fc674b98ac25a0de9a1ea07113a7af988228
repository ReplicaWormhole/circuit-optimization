"""Bounded CNOT-power homotopy from exact14 to 13 CNOTs.

Replace one exact14 CX(c,t) by H_t CP(pi*alpha) H_t. At alpha=1 this is
exactly CX(c,t); at alpha=0 it is identity. At fixed descending alpha
stages, refit arbitrary SU(2) layers between every entangler using the
label-free squared off-diagonal Frobenius loss of U V4 U-dagger /16.
Only the alpha=0 endpoint has a 13-CNOT gate list. Failed paths give no
lower bound, and a numerical near-hit would still need exactification.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import embedded_one_qubit, evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
POSITIONS = (3, 6, 12)
STAGES = (0.875, 0.75, 0.5, 0.25, 0.0)
SIGMAS = (0.02, 0.15)
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64


def split_exact14():
    source = json.loads(BASE.read_text())
    chunks, cxs = [[]], []
    for gate in source["gates"]:
        if gate["gate"] == "cx":
            cxs.append(dict(gate))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if len(cxs) != 14 or len(chunks) != 15:
        raise ValueError("expected exact14 source with 14 CX and 15 local chunks")
    fixed = []
    for chunk in chunks:
        u = np.eye(16, dtype=complex)
        for gate in chunk:
            u = fixed_matrix(gate) @ u
        fixed.append(torch.as_tensor(u, dtype=CDTYPE))
    return chunks, cxs, fixed


def powered_cx(gate, alpha):
    c, t = gate["control"], gate["target"]
    h = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
    ht = embedded_one_qubit(h, t, 4)
    phases = [np.exp(1j * math.pi * alpha)
              if (basis & (8 >> c)) and (basis & (8 >> t)) else 1
              for basis in range(16)]
    return torch.as_tensor(ht @ np.diag(phases) @ ht, dtype=CDTYPE)


def fit_stage(x0, alpha, removed, cxs, fixed, maxiter):
    cxm = [powered_cx(g, alpha if j == removed else 1.0)
           for j, g in enumerate(cxs)]
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)

    def objective(flat):
        angles = torch.tensor(flat.reshape(15, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(15):
            u = fixed[slot] @ u
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
            if slot < 14:
                u = cxm[slot] @ u
        d = u @ v @ u.conj().T
        off = d - torch.diag(torch.diagonal(d))
        loss = (off.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(x0.ravel())[0]
    opt = minimize(objective, x0.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    return opt.x.reshape(15, 4, 3), {
        "alpha": alpha, "initial_loss": initial, "loss": float(opt.fun),
        "nit": int(opt.nit), "nfev": int(opt.nfev),
        "optimizer_success": bool(opt.success),
        "optimizer_message": str(opt.message),
    }


def serialize_endpoint(chunks, cxs, angles, removed):
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < 14 and slot != removed:
            gates.append(cxs[slot])
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=14200)
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--prefix", default="search13_cx_power_homotopy")
    args = parser.parse_args()
    chunks, cxs, fixed = split_exact14()
    matrix_error = max(float(torch.max(torch.abs(powered_cx(g, 1.0)
                       - torch.as_tensor(fixed_matrix(g), dtype=CDTYPE))))
                       for g in cxs)
    if matrix_error > 1e-12:
        raise AssertionError("alpha=1 gate does not equal CX")
    rows = []
    for pos_index, removed in enumerate(POSITIONS):
        for start, sigma in enumerate(SIGMAS):
            seed = args.seed + 2 * pos_index + start
            rng = np.random.default_rng(seed)
            angles = rng.normal(0, sigma, (15, 4, 3))
            stages = []
            for alpha in STAGES:
                angles, row = fit_stage(angles, alpha, removed, cxs,
                                        fixed, args.maxiter)
                stages.append(row)
                print(json.dumps({"removed": removed, "seed": seed,
                                  "stage": row}, sort_keys=True), flush=True)
            candidate = serialize_endpoint(chunks, cxs, angles, removed)
            path = ROOT / f"{args.prefix}_p{removed}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            check = evaluate(candidate)
            record = {"removed": removed, "seed": seed, "sigma": sigma,
                      "maxiter_per_stage": args.maxiter, "stages": stages,
                      "candidate": path.name,
                      "endpoint_off_diagonal_error": check["off_diagonal_error"],
                      "endpoint_eigenvalue_root_error": check[
                          "eigenvalue_root_error"],
                      "endpoint_valid_diagonalizer": check["valid_diagonalizer"]}
            rows.append(record)
            print(json.dumps({"endpoint": record}, sort_keys=True), flush=True)
    result = {"base": BASE.name, "positions": POSITIONS, "stages": STAGES,
              "sigmas": SIGMAS, "alpha1_cx_matrix_error": matrix_error,
              "rows": rows, "best": min(rows, key=lambda r:
              r["stages"][-1]["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name,
                      "best_candidate": result["best"]["candidate"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
