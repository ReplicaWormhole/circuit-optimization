"""Bounded fixed-unitary seven-CNOT fits for a same-eigenvalue row swap.

Rows 8 and 11 of the exact14 diagonalizer both have V4 eigenvalue +1.
Swapping them preserves diagonalization and lowers one balanced cut rank of
the full tail from 16 to 12. Two nearest rank-compatible seven-CNOT graphs
are fitted with arbitrary SU(2) layers around all CNOTs. This is a numeric
screen; any low-loss candidate needs independent exactification.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
SCHEDULES = (
    ((0, 2), (0, 2), (1, 3), (0, 3), (0, 1), (0, 1), (2, 3)),
    ((0, 2), (0, 2), (1, 3), (0, 3), (0, 1), (2, 3), (2, 3)),
)
SEEDS = (2201, 2202)
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64


def target_and_prefix():
    gates = json.loads(BASE.read_text())["gates"]
    prefix, tail = gates[:13], gates[13:]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    assert sum(g["gate"] == "cx" for g in tail) == 8
    target = np.eye(16, dtype=complex)
    for gate in tail:
        target = fixed_matrix(gate) @ target
    target[[8, 11]] = target[[11, 8]]
    return prefix, target


def fit(edges, seed, maxiter, target):
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0, 0.4, (8, 4, 3))
    target_t = torch.as_tensor(target, dtype=CDTYPE)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": c,
                                         "target": t}), dtype=CDTYPE)
           for c, t in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(8):
            for qubit in range(4):
                u = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ u
            if slot < 7:
                u = cxs[slot] @ u
        overlap = torch.trace(target_t.conj().T @ u)
        loss = 1 - overlap.abs().square().real / 256
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(x0.ravel())[0]
    opt = minimize(objective, x0.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    angles = opt.x.reshape(8, 4, 3)
    gates = []
    for slot in range(8):
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 7:
            control, target_wire = edges[slot]
            gates.append({"gate": "cx", "control": control,
                          "target": target_wire})
    return gates, {"seed": seed, "maxiter": maxiter,
                   "initial_loss": initial, "loss": float(opt.fun),
                   "nit": int(opt.nit), "nfev": int(opt.nfev),
                   "optimizer_success": bool(opt.success),
                   "message": str(opt.message)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=250)
    args = parser.parse_args()
    prefix, target = target_and_prefix()
    rows = []
    for schedule_index, edges in enumerate(SCHEDULES):
        for seed in SEEDS:
            tail, record = fit(edges, seed, args.maxiter, target)
            candidate = {"n": 4, "gates": list(prefix) + tail}
            path = ROOT / f"search13_eigenrow_gauge_seven_fit_top{schedule_index}_seed{seed}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            check = evaluate(candidate)
            record.update(schedule_index=schedule_index, edges=edges,
                          candidate=path.name,
                          off_diagonal_error=check["off_diagonal_error"],
                          valid_diagonalizer=check["valid_diagonalizer"])
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    result = {"gauge": {"kind": "same-eigenvalue-row-swap", "rows": [8, 11]},
              "schedules": SCHEDULES, "seeds": SEEDS,
              "maxiter": args.maxiter, "rows": rows,
              "best": min(rows, key=lambda row: row["loss"])}
    out = ROOT / "search13_eigenrow_gauge_seven_fit_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": out.name, "best": result["best"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
