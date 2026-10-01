"""Bounded 5+8 boundary-bridge synthesis of the exact rank-eight gauge.

The sixth CNOT of exact14's prefix is moved into the middle of the tail and
one tail CNOT is removed.  All one-qubit gates are then optimized.  The first
stage targets the exact gauged unitary up to global phase; the second uses the
label-free cycle condition.  Numerical fits are not exact certificates.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, one_qubit_matrix, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate
from search13_joint_gauge_five_endpoint import local_to_axis


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_rank8_bridge_result.json"
EDGES = ((3, 1), (1, 2), (3, 2), (0, 2), (1, 2),
         (0, 2), (0, 2), (1, 3), (1, 3), (3, 2),
         (0, 1), (0, 1), (2, 3))
torch.set_num_threads(1)


def exact_gauge():
    h = 1 / np.sqrt(2)
    g = np.zeros((16, 16), dtype=complex)
    entries = {
        0: {0: 1}, 1: {1: 1j*h, 14: h}, 2: {2: 1}, 3: {13: 1},
        4: {6: -h, 9: h}, 5: {5: h, 10: h}, 6: {4: 1}, 7: {7: 1},
        8: {15: 1}, 9: {6: h, 9: h}, 10: {5: -h, 10: h},
        11: {11: 1}, 12: {12: 1}, 13: {3: 1},
        14: {1: 1j*h, 14: -h}, 15: {8: 1},
    }
    for row, columns in entries.items():
        for col, value in columns.items():
            g[row, col] = value
    assert np.max(np.abs(g @ g.conj().T - np.eye(16))) < 1e-14
    return g


def baseline_and_warm():
    gates = json.loads(SOURCE.read_text())["gates"]
    u = np.eye(16, dtype=complex)
    for gate in gates:
        u = fixed_matrix(gate) @ u
    chunks = [[]]
    for gate in gates:
        if gate["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(gate)
    # Delete the last CNOT in the 14-CNOT source, then combine its adjacent
    # local layers.  This supplies a reproducible 13-CNOT warm parameter set.
    assert len(chunks) == 15
    chunks[13].extend(chunks[14])
    chunks = chunks[:14]
    layers = []
    for chunk in chunks:
        mats = [np.eye(2, dtype=complex) for _ in range(4)]
        for gate in chunk:
            q = gate["qubit"]
            mats[q] = one_qubit_matrix(gate) @ mats[q]
        layers.append([local_to_axis(mat) for mat in mats])
    return u, np.asarray(layers)


def matrix(angles, cx_matrices):
    u = torch.eye(16, dtype=torch.complex128)
    for slot in range(14):
        for wire in range(4):
            u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
        if slot < 13:
            u = cx_matrices[slot] @ u
    return u


def optimize(initial, target, shift, cxs, mode, maxiter):
    def objective(flat):
        angles = torch.tensor(flat.reshape(14, 4, 3), dtype=torch.float64,
                              requires_grad=True)
        u = matrix(angles, cxs)
        if mode == "process":
            overlap = torch.trace(target.conj().T @ u)
            loss = 1 - overlap.abs().square() / 256
        else:
            a = u @ shift @ u.conj().T
            off = a - torch.diag(torch.diagonal(a))
            loss = off.abs().square().sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()
    start_loss = objective(initial.ravel())[0]
    fit = minimize(objective, initial.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    return fit.x.reshape(14, 4, 3), {
        "start_loss": start_loss, "loss": float(fit.fun),
        "iterations": int(fit.nit), "success": bool(fit.success),
        "message": str(fit.message),
    }


def candidate(angles):
    gates = []
    for slot in range(14):
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 13:
            a, b = EDGES[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=26613)
    parser.add_argument("--process-maxiter", type=int, default=150)
    parser.add_argument("--cycle-maxiter", type=int, default=200)
    args = parser.parse_args()
    exact_u, warm = baseline_and_warm()
    g = exact_gauge()
    target = torch.as_tensor(g @ exact_u, dtype=torch.complex128)
    shift = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                        "target": b}), dtype=torch.complex128)
           for a, b in EDGES]
    rng = np.random.default_rng(args.seed)
    rows = []
    for name, init in (("warm", warm + rng.normal(0, 0.03, warm.shape)),
                       ("random", rng.normal(0, 0.4, warm.shape))):
        angles, process = optimize(init, target, shift, cxs,
                                   "process", args.process_maxiter)
        angles, cycle = optimize(angles, target, shift, cxs,
                                 "cycle", args.cycle_maxiter)
        path = ROOT / f"search13_rank8_bridge_{name}.json"
        circuit = candidate(angles)
        path.write_text(json.dumps(circuit, indent=2) + "\n")
        check = evaluate(circuit)
        row = {"start": name, "process": process, "cycle": cycle,
               "candidate": path.name,
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    OUT.write_text(json.dumps({
        "source": SOURCE.name, "gauge_source":
            "search13_fulltail_invariant_exact_witness.py",
        "edges": EDGES, "seed": args.seed,
        "process_maxiter": args.process_maxiter,
        "cycle_maxiter": args.cycle_maxiter, "rows": rows,
        "scope": "two bounded numerical fits of one boundary-bridge topology",
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
