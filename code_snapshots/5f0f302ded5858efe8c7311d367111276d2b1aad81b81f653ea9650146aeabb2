"""Bounded nonmonomial five-CNOT prefix and eight-CNOT tail homotopy.

Use an endpoint-guided path not fitted in run 249.  Seed transverse rotations
inside the prefix, then optimize every local SU(2) layer.  Stages target the
exact14 output basis by process overlap, its eigenvalue assignment by an
eigen-residual, and finally the label-free cycle condition.  Numerical fits
cannot prove impossibility.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate
from search13_endpoint_guided_fit import phase_informed_angles, tail_chunks_and_edges
from search13_joint_gauge_five_endpoint import endpoint


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
SCREEN = ROOT / "search13_endpoint_guided_screen_result.json"
OUT = ROOT / "search13_nonmonomial_homotopy_result.json"
SELECTED_ORIENTATION = "2->0"
SELECTED_RANK = 1  # Run 249 used 0 and 3, not this path.
torch.set_num_threads(1)


def circuit_matrix(gates):
    u = np.eye(16, dtype=complex)
    for gate in gates:
        u = fixed_matrix(gate) @ u
    return u


def unitary(angles, cxs):
    u = torch.eye(16, dtype=torch.complex128)
    for slot in range(14):
        for wire in range(4):
            u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
        if slot < 13:
            u = cxs[slot] @ u
    return u


def fit_stage(initial, cxs, target, shift, diagonal, mode, maxiter):
    def objective(flat):
        angles = torch.tensor(flat.reshape(14, 4, 3), dtype=torch.float64,
                              requires_grad=True)
        u = unitary(angles, cxs)
        if mode == "process":
            overlap = torch.trace(target.conj().T @ u)
            loss = 1 - overlap.abs().square() / 256
        elif mode == "fixed-label":
            residual = u @ shift - diagonal @ u
            loss = residual.abs().square().sum().real / 16
        else:
            transformed = u @ shift @ u.conj().T
            off = transformed - torch.diag(torch.diagonal(transformed))
            loss = off.abs().square().sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    start_loss = objective(initial.ravel())[0]
    opt = minimize(objective, initial.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    return opt.x.reshape(14, 4, 3), {
        "mode": mode, "start_loss": start_loss, "loss": float(opt.fun),
        "iterations": int(opt.nit), "success": bool(opt.success),
        "message": str(opt.message),
    }


def candidate(edges, angles):
    gates = []
    for slot in range(14):
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 13:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return {"n": 4, "gates": gates}


def nonmonomial_score(gates):
    prefix = []
    count = 0
    for gate in gates:
        prefix.append(gate)
        if gate["gate"] == "cx":
            count += 1
            if count == 5:
                break
    u = circuit_matrix(prefix)
    probabilities = np.abs(u)**2
    return float(np.max(1 - np.max(probabilities, axis=0)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=27913)
    parser.add_argument("--process-maxiter", type=int, default=120)
    parser.add_argument("--fixed-maxiter", type=int, default=180)
    parser.add_argument("--cycle-maxiter", type=int, default=200)
    args = parser.parse_args()
    screen = json.loads(SCREEN.read_text())
    selected = screen["selected"][SELECTED_ORIENTATION][SELECTED_RANK]
    prefix = tuple(tuple(edge) for edge in selected["path"])
    assert endpoint(prefix) == tuple(selected["endpoint"])
    assert selected["deletion_distance"] >= 2
    tail_chunks, tail_edges = tail_chunks_and_edges()
    edges = prefix + tail_edges
    seeded, applied = phase_informed_angles(prefix, tail_chunks)
    assert len(edges) == 13 and len(applied) == 2
    source = json.loads(BASE.read_text())
    target_np = circuit_matrix(source["gates"])
    shift_np = right_shift(4)
    d = target_np @ shift_np @ target_np.conj().T
    assert np.max(np.abs(d - np.diag(np.diag(d)))) < 1e-12
    target = torch.as_tensor(target_np, dtype=torch.complex128)
    shift = torch.as_tensor(shift_np, dtype=torch.complex128)
    diagonal = torch.diag(torch.as_tensor(np.diag(d), dtype=torch.complex128))
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                        "target": b}), dtype=torch.complex128)
           for a, b in edges]
    rng = np.random.default_rng(args.seed)
    rows = []
    for sign in (1, -1):
        initial = seeded.copy()
        initial[0, 2, 0] += sign * 0.55
        initial[2, 1, 1] += sign * 0.55
        initial += rng.normal(0, 0.025, initial.shape)
        initial_nonmonomial = nonmonomial_score(candidate(edges, initial)["gates"])
        angles = initial
        stages = []
        for mode, cap in (("process", args.process_maxiter),
                          ("fixed-label", args.fixed_maxiter),
                          ("cycle", args.cycle_maxiter)):
            angles, row = fit_stage(angles, cxs, target, shift,
                                    diagonal, mode, cap)
            stages.append(row)
        circuit = candidate(edges, angles)
        name = "plus" if sign > 0 else "minus"
        path = ROOT / f"search13_nonmonomial_homotopy_{name}.json"
        path.write_text(json.dumps(circuit, indent=2) + "\n")
        checked = evaluate(circuit)
        row = {"start": name, "initial_nonmonomial_score": initial_nonmonomial,
               "final_nonmonomial_score": nonmonomial_score(circuit["gates"]),
               "stages": stages, "candidate": path.name,
               "off_diagonal_error": checked["off_diagonal_error"],
               "valid_diagonalizer": checked["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    OUT.write_text(json.dumps({
        "source": BASE.name, "screen": SCREEN.name,
        "selected_orientation": SELECTED_ORIENTATION,
        "selected_rank": SELECTED_RANK, "prefix": prefix,
        "endpoint": selected["endpoint"], "tail": tail_edges,
        "visited_exact_phase_masks": applied, "seed": args.seed,
        "process_maxiter": args.process_maxiter,
        "fixed_maxiter": args.fixed_maxiter,
        "cycle_maxiter": args.cycle_maxiter, "rows": rows,
        "scope": "one previously unfitted 5+8 topology, two transverse seeds, "
                 "three-stage numerical homotopy",
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
