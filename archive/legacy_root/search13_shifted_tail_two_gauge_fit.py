"""Bounded seven-CNOT synthesis fit of a gauged exact14 full tail.

The target is G*T, with T the eight-CNOT exact14 tail after its six-CNOT
prefix and G=exp(-i*pi/8*Z_5) exp(i*pi/8*Z_6). Two distinct seven-CNOT
star-plus-chord schedules pass the exact cut-rank graph screen (run 209).
Fit arbitrary SU(2) layers before, between, and after CNOTs by maximizing
phase-insensitive process overlap. A numerical fit is not an exact circuit.
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
    ((0, 2), (0, 2), (0, 3), (0, 3), (0, 1), (0, 1), (1, 2)),
    ((0, 1), (0, 2), (0, 3), (0, 1), (0, 2), (0, 3), (1, 2)),
)
MASKS = (5, 6)
ANGLE_NUMERATORS = (-2, 2)
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
    phase = np.array([np.exp(1j * np.pi / 16 * sum(
        k * (1 if (basis & mask).bit_count() % 2 == 0 else -1)
        for mask, k in zip(MASKS, ANGLE_NUMERATORS)))
        for basis in range(16)])
    return prefix, phase[:, None] * target


def fit(edges, seed, sigma, maxiter, target):
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0, sigma, (8, 4, 3))
    target_t = torch.as_tensor(target, dtype=CDTYPE)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": c,
                                          "target": t}), dtype=CDTYPE)
           for c, t in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(8):
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
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
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < 7:
            c, t = edges[slot]
            gates.append({"gate": "cx", "control": c, "target": t})
    return gates, {"seed": seed, "sigma": sigma, "maxiter": maxiter,
                   "initial_loss": initial, "loss": float(opt.fun),
                   "nit": int(opt.nit), "nfev": int(opt.nfev),
                   "optimizer_success": bool(opt.success),
                   "optimizer_message": str(opt.message)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=13900)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--prefix", default="search13_shifted_tail_two_gauge_fit")
    args = parser.parse_args()
    prefix, target = target_and_prefix()
    rows = []
    for schedule_index, edges in enumerate(SCHEDULES):
        for start, sigma in enumerate((0.3, 1.0)):
            seed = args.seed + 2 * schedule_index + start
            tail_gates, record = fit(edges, seed, sigma, args.maxiter, target)
            candidate = {"n": 4, "gates": list(prefix) + tail_gates}
            path = ROOT / f"{args.prefix}_top{schedule_index}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            check = evaluate(candidate)
            record.update(schedule_index=schedule_index,
                          edges=edges, candidate=path.name,
                          off_diagonal_error=check["off_diagonal_error"],
                          eigenvalue_root_error=check["eigenvalue_root_error"],
                          valid_diagonalizer=check["valid_diagonalizer"])
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    summary = {"base": BASE.name, "boundary_gate_start": 13,
               "masks": MASKS, "angle_numerators": ANGLE_NUMERATORS,
               "schedules": SCHEDULES, "rows": rows,
               "best": min(rows, key=lambda r: r["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "best": summary["best"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
