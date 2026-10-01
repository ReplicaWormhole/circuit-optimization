"""Bounded 13-CNOT matchgate deletion plus neighboring-edge reroute.

The six prescribed topologies delete the first CNOT of one exact matchgate
block and reroute its surviving CNOT to an edge of an adjacent block.  Each
topology gets two independent SU(2)-layer initializations.  The loss is the
label-free squared Frobenius norm of the off-diagonal part of U V4 U† / 16.
This script is a numerical search, not an exact certificate.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
CASES = [
    (6, 7, (1, 3)),
    (8, 9, (0, 2)),
    (8, 9, (0, 1)),
    (10, 11, (1, 3)),
    (10, 11, (2, 3)),
    (12, 13, (0, 1)),
]

torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64


def topology(deleted, rerouted, edge):
    source = json.loads(BASE.read_text())
    chunks, cxs = [[]], []
    old_to_new = {}
    index = -1
    for gate in source["gates"]:
        if gate["gate"] == "cx":
            index += 1
            if index == deleted:
                continue
            old_to_new[index] = len(cxs)
            cxs.append(dict(gate))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if index != 13 or len(cxs) != 13:
        raise ValueError("expected an exact 14-CNOT source")
    cxs[old_to_new[rerouted]].update(control=edge[0], target=edge[1])
    fixed = []
    for chunk in chunks:
        matrix = np.eye(16, dtype=complex)
        for gate in chunk:
            matrix = fixed_matrix(gate) @ matrix
        fixed.append(torch.as_tensor(matrix, dtype=CDTYPE))
    cxm = [torch.as_tensor(fixed_matrix(gate), dtype=CDTYPE) for gate in cxs]
    return chunks, cxs, fixed, cxm


def optimize(deleted, rerouted, edge, seed, sigma, maxiter):
    chunks, cxs, fixed, cxm = topology(deleted, rerouted, edge)
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0, sigma, (14, 4, 3))

    def objective(flat):
        angles = torch.tensor(flat.reshape(14, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(14):
            u = fixed[slot] @ u
            for qubit in range(4):
                u = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ u
            if slot < 13:
                u = cxm[slot] @ u
        diagonalized = u @ v @ u.conj().T
        off = diagonalized - torch.diag(torch.diagonal(diagonalized))
        loss = (off.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    opt = minimize(objective, x0.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    angles = opt.x.reshape(14, 4, 3)
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 13:
            gates.append(cxs[slot])
    candidate = {"n": 4, "gates": gates}
    check = evaluate(candidate)
    return candidate, {
        "deleted": deleted, "rerouted": rerouted, "edge": list(edge),
        "seed": seed, "sigma": sigma, "maxiter": maxiter,
        "nit": int(opt.nit), "nfev": int(opt.nfev),
        "optimizer_success": bool(opt.success),
        "optimizer_message": str(opt.message),
        "loss": float(opt.fun),
        "off_diagonal_error": check["off_diagonal_error"],
        "eigenvalue_root_error": check["eigenvalue_root_error"],
        "valid_diagonalizer": check["valid_diagonalizer"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--seed", type=int, default=13150)
    parser.add_argument("--prefix", type=str,
                        default="search13_reroute_neighbor")
    args = parser.parse_args()
    results = []
    for case_index, (deleted, rerouted, edge) in enumerate(CASES):
        for start in range(2):
            seed = args.seed + 2 * case_index + start
            sigma = 0.04 if start == 0 else 0.3
            candidate, record = optimize(deleted, rerouted, edge, seed,
                                         sigma, args.maxiter)
            path = ROOT / f"{args.prefix}_d{deleted}_r{rerouted}_{edge[0]}{edge[1]}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            record["candidate"] = path.name
            results.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    summary = {"base": BASE.name, "results": results,
               "best": min(results, key=lambda r: r["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "best": summary["best"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
