"""Bounded numerical search for a 16-CNOT Bell-coordinate V4 diagonalizer.

Chronological circuit: Bell-coordinate inverse on (0,2),(1,3), compute the
two-bit label XOR, then optimize 12 entanglers and local SU(2) layers.
This is deliberately distinct from deleting a gate from the 18-CNOT circuit.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch

from check_circuit import cnot, right_shift


torch.set_num_threads(1)
DTYPE = torch.complex128
I2 = torch.eye(2, dtype=DTYPE)
PAULI_Z = torch.diag(torch.tensor([1.0, -1.0], dtype=DTYPE))
PAULI_Y = torch.tensor([[0, -1j], [1j, 0]], dtype=DTYPE)
H = torch.tensor([[1, 1], [1, -1]], dtype=DTYPE) / math.sqrt(2)


def embed(a, q):
    out = torch.tensor([[1.0]], dtype=DTYPE)
    for j in range(4):
        out = torch.kron(out, a if j == q else I2)
    return out


H0, H1 = embed(H, 0), embed(H, 1)
CACHE = {(c, t): torch.tensor(cnot(c, t, 4), dtype=DTYPE)
         for c in range(4) for t in range(4) if c != t}
V = torch.tensor(right_shift(4), dtype=DTYPE)
PREFIX = [(0, 2), (1, 3), (0, 1), (2, 3)]
PAIRS = [(0, 2), (0, 1), (2, 3), (0, 3), (2, 1), (1, 3)]


def topology(seed):
    # Paired entanglers initialize to identity, then local gates can deform them.
    rng = np.random.default_rng(seed)
    sweep = PAIRS.copy()
    rng.shuffle(sweep)
    return [pair for pair0 in sweep for pair in (pair0, pair0)]


def local_matrix(angles):
    mats = []
    for q in range(4):
        a, b, c = angles[q]
        za = torch.cos(a / 2) * I2 - 1j * torch.sin(a / 2) * PAULI_Z
        yb = torch.cos(b / 2) * I2 - 1j * torch.sin(b / 2) * PAULI_Y
        zc = torch.cos(c / 2) * I2 - 1j * torch.sin(c / 2) * PAULI_Z
        mats.append(zc @ yb @ za)
    out = mats[0]
    for m in mats[1:]:
        out = torch.kron(out, m)
    return out


def unitary(angles, suffix):
    u = torch.eye(16, dtype=DTYPE)
    # B†: CX then H on the first qubit of each Bell pair.
    u = H0 @ CACHE[PREFIX[0]] @ u
    u = H1 @ CACHE[PREFIX[1]] @ u
    u = CACHE[PREFIX[2]] @ u
    u = CACHE[PREFIX[3]] @ u
    for i, pair in enumerate(suffix):
        u = CACHE[pair] @ local_matrix(angles[i]) @ u
    return local_matrix(angles[-1]) @ u


def loss(angles, suffix):
    u = unitary(angles, suffix)
    m = u @ V @ u.conj().T
    off = m - torch.diag(torch.diag(m))
    return torch.sum(torch.abs(off) ** 2) / 16


def gates(angles, suffix):
    out = [dict(gate="cx", control=0, target=2), dict(gate="h", qubit=0),
           dict(gate="cx", control=1, target=3), dict(gate="h", qubit=1),
           dict(gate="cx", control=0, target=1),
           dict(gate="cx", control=2, target=3)]
    for i in range(len(suffix) + 1):
        for q in range(4):
            for name, theta in zip(("rz", "ry", "rz"), angles[i, q]):
                out.append(dict(gate=name, qubit=q, theta=float(theta)))
        if i < len(suffix):
            c, t = suffix[i]
            out.append(dict(gate="cx", control=c, target=t))
    return {"n": 4, "gates": out}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--lr", type=float, default=0.06)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    suffix = topology(args.seed)
    initial = 0.025 * torch.randn(13, 4, 3, dtype=torch.float64)
    # The base circuit H_a XOR B† diagonalizes unsigned pair-register SWAP.
    initial[-1, 0, :] += torch.tensor([math.pi, math.pi / 2, 0])
    initial[-1, 2, :] += torch.tensor([math.pi, math.pi / 2, 0])
    angles = torch.nn.Parameter(initial)
    optimizer = torch.optim.Adam([angles], lr=args.lr)
    best = float("inf")
    best_angles = None
    for step in range(args.steps):
        optimizer.zero_grad()
        objective = loss(angles, suffix)
        objective.backward()
        optimizer.step()
        value = float(objective.detach())
        if value < best:
            best, best_angles = value, angles.detach().clone().numpy()
        if step % 100 == 0:
            print(json.dumps({"step": step, "loss": value, "best": best}), flush=True)
    candidate = gates(best_angles, suffix)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    print(json.dumps({"best_frobenius_squared_over_16": best,
                      "suffix_topology": suffix, "candidate": str(args.output)}))


if __name__ == "__main__":
    main()
