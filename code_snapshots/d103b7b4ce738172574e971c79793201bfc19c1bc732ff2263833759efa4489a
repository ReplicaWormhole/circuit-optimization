"""Bounded 13-CNOT deletion/reoptimization from the exact 14-CX candidate."""
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate

torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
DELETIONS = [6, 7, 8, 9, 10, 11, 12, 13]


def run(deleted, seed, maxiter, out):
    source = json.loads(BASE.read_text())
    chunks = [[]]
    cxs = []
    ci = -1
    for gate in source["gates"]:
        if gate["gate"] == "cx":
            ci += 1
            if ci == deleted:
                continue
            cxs.append(gate)
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if ci != 13 or len(cxs) != 13:
        raise ValueError("source must have 14 CNOTs")
    fixed = []
    for chunk in chunks:
        u = np.eye(16, dtype=complex)
        for gate in chunk:
            u = fixed_matrix(gate) @ u
        fixed.append(torch.as_tensor(u, dtype=CDTYPE))
    cxm = [torch.as_tensor(fixed_matrix(g), dtype=CDTYPE) for g in cxs]
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0, 0.04, (14, 4, 3))

    def objective(flat):
        angles = torch.tensor(flat.reshape(14, 4, 3), dtype=RDTYPE, requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(14):
            u = fixed[slot] @ u
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
            if slot < 13:
                u = cxm[slot] @ u
        transformed = u @ v @ u.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel()

    opt = minimize(objective, x0.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-15, "gtol": 1e-11,
                            "maxls": 30})
    angles = opt.x.reshape(14, 4, 3)
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < 13:
            gates.append(cxs[slot])
    candidate = {"n": 4, "gates": gates}
    out.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    print(json.dumps({"deleted_cx": deleted, "seed": seed,
                      "maxiter": maxiter, "nit": opt.nit, "nfev": opt.nfev,
                      "objective": opt.fun, "success": bool(opt.success),
                      "message": str(opt.message),
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid_diagonalizer": check["valid_diagonalizer"],
                      "candidate": str(out)}, sort_keys=True))


if __name__ == "__main__":
    for index in DELETIONS:
        run(index, seed=1300 + index, maxiter=250,
            out=ROOT / f"search13_matchgate_delete_cx{index}.json")
