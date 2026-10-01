"""Bounded joint seven-CNOT tail fit from the exact 14-CNOT circuit.

Fix the six-CNOT parity prefix, remove one CNOT from the first matchgate
layer, reroute another CNOT, and optimize every local SU(2) frame from the
tail input through its output. The label-free loss permits eigenspace gauge.
This is numerical search, not an optimality test.
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


def prepare(deleted, rerouted, edge):
    source = json.loads(BASE.read_text())
    prefix, tail = source["gates"][:13], source["gates"][13:]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    chunks, edges = [[]], []
    for gate in tail:
        if gate["gate"] == "cx":
            edges.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(edges) == 8
    if deleted not in range(4) or rerouted not in range(4) or deleted == rerouted:
        raise ValueError("delete and reroute distinct first-layer CNOT indices 0..3")
    edges[rerouted] = tuple(edge)
    del edges[deleted]
    chunks[deleted].extend(chunks.pop(deleted + 1))
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
    return prefix, chunks, edges, mats, cmats, torch.as_tensor(prefix_u, dtype=CDTYPE)


def run(args):
    prefix, chunks, edges, mats, cmats, prefix_u = prepare(
        args.delete, args.reroute, args.edge)
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    rng = np.random.default_rng(args.seed)
    x0 = rng.normal(0, args.sigma, (8, 4, 3))

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = prefix_u
        for slot in range(8):
            u = mats[slot] @ u
            for q in range(4):
                u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
            if slot < 7:
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
    angles = opt.x.reshape(8, 4, 3)
    gates = list(prefix)
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < 7:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    candidate = {"n": 4, "gates": gates}
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    print(json.dumps({"delete": args.delete, "reroute": args.reroute,
                      "edge": args.edge, "seed": args.seed,
                      "maxiter": args.maxiter, "initial_loss": initial,
                      "final_loss": opt.fun, "nit": opt.nit,
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid_diagonalizer": check["valid_diagonalizer"],
                      "topology": edges, "candidate": str(args.output)},
                     sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--delete", type=int, required=True)
    parser.add_argument("--reroute", type=int, required=True)
    parser.add_argument("--edge", type=int, nargs=2, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--sigma", type=float, default=0.04)
    parser.add_argument("--maxiter", type=int, default=400)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
