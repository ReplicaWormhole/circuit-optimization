"""Fit the exact-14 four-matchgate topology after a new five-CNOT parity prefix.

Unlike earlier five-prefix scans, the suffix topology is the eight-CNOT
four-matchgate tail (02,13,01,23). All nine local tail frames vary and the
label-free loss allows degenerate eigenspace mixing. The tested five-step
parity tours use a different diagonal gauge/endpoint from the exact 14-CNOT
prefix. A failed fit excludes only the stated initializations and topology.
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
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64

# The cancelled triple has zero phase after shifting all triple coefficients
# by an invariant diagonal gauge. The tours expose the other required masks.
TRIPLES = (11, 13, 14, 7)
COEFFICIENTS = (1, 2, 3, 0)
SINGLES = {1: "5*pi/8", 2: "3*pi/4", 3: "3*pi/8"}
TOURS = {
    "g0_near": (0, ((0, 2), (2, 1), (0, 1), (3, 1), (2, 1))),
    "g2_near_a": (2, ((2, 1), (0, 2), (3, 1), (2, 1), (2, 3))),
    "g2_near_b": (2, ((1, 2), (0, 1), (3, 2), (1, 2), (3, 1))),
}


def prefix_gates(name):
    cancelled, tour = TOURS[name]
    shift = -COEFFICIENTS[cancelled]
    phases = {mask: f"{coefficient + shift}*pi/8"
              for mask, coefficient in zip(TRIPLES, COEFFICIENTS)
              if coefficient + shift}
    phases[6] = "-pi/2"
    gates = [{"gate": "rz", "qubit": q, "theta": angle}
             for q, angle in SINGLES.items()]
    rows, seen = [8, 4, 2, 1], set()
    for c, t in tour:
        gates.append({"gate": "cx", "control": c, "target": t})
        rows[t] ^= rows[c]
        mask = rows[t]
        if mask in phases and mask not in seen:
            gates.append({"gate": "rz", "qubit": t, "theta": phases[mask]})
            seen.add(mask)
    assert seen == set(phases)
    return gates, tuple(rows), phases


def prepare(name):
    prefix, endpoint, phases = prefix_gates(name)
    source = json.loads(BASE.read_text())
    tail = source["gates"][13:]
    chunks, edges = [[]], []
    for gate in tail:
        if gate["gate"] == "cx":
            edges.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert edges == [(0, 2)] * 2 + [(1, 3)] * 2 + [(0, 1)] * 2 + [(2, 3)] * 2
    mats = []
    for chunk in chunks:
        u = np.eye(16, dtype=complex)
        for gate in chunk:
            u = fixed_matrix(gate) @ u
        mats.append(torch.as_tensor(u, dtype=CDTYPE))
    cmats = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                           "target": b}), dtype=CDTYPE)
             for a, b in edges]
    prefix_u = np.eye(16, dtype=complex)
    for gate in prefix:
        prefix_u = fixed_matrix(gate) @ prefix_u
    return prefix, endpoint, phases, chunks, edges, mats, cmats, torch.as_tensor(prefix_u, dtype=CDTYPE)


def run(args):
    prefix, endpoint, phases, chunks, edges, mats, cmats, prefix_u = prepare(args.tour)
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    rng = np.random.default_rng(args.seed)
    x0 = rng.normal(0, args.sigma, (9, 4, 3))

    def objective(flat):
        angles = torch.tensor(flat.reshape(9, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = prefix_u
        for slot in range(9):
            u = mats[slot] @ u
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
            if slot < 8:
                u = cmats[slot] @ u
        transformed = u @ v @ u.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(x0.ravel())[0]
    opt = minimize(objective, x0.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": args.maxiter, "ftol": 1e-15,
                            "gtol": 1e-11, "maxls": 30})
    angles = opt.x.reshape(9, 4, 3)
    gates = list(prefix)
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < 8:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    candidate = {"n": 4, "gates": gates}
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    print(json.dumps({"tour": args.tour, "prefix_endpoint": endpoint,
                      "phase_masks": phases, "seed": args.seed,
                      "maxiter": args.maxiter, "initial_loss": initial,
                      "final_loss": opt.fun, "nit": opt.nit,
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid_diagonalizer": check["valid_diagonalizer"],
                      "candidate": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tour", choices=TOURS, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--sigma", type=float, default=0.04)
    parser.add_argument("--maxiter", type=int, default=600)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
