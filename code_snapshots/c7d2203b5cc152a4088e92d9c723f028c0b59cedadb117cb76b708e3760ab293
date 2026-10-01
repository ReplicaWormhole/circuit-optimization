"""Optimize four-CNOT connected tails after the exact ten-CNOT prefix.

The loss is the full off-diagonal Frobenius norm of U V4 U†, so the search
permits any output eigenvalue order and any basis in degenerate eigenspaces.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from check_circuit import (cnot, embedded_one_qubit, one_qubit_matrix,
                           right_shift)
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)
DTYPE = torch.complex128
I2 = torch.eye(2, dtype=DTYPE)
Z = torch.diag(torch.tensor([1.0, -1.0], dtype=DTYPE))
Y = torch.tensor([[0, -1j], [1j, 0]], dtype=DTYPE)
V = torch.tensor(right_shift(4), dtype=DTYPE)
TOPOLOGIES = (
    ((2, 1), (0, 1), (2, 3), (0, 3)),
    ((1, 3), (2, 1), (0, 2), (0, 1)),
    ((0, 1), (1, 2), (2, 3), (0, 1)),
    ((2, 1), (2, 3), (2, 0), (0, 1)),
)


def prefix_data():
    candidate = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    split = 13 + len(exact_block()) + len(reverse_block())
    gates = candidate["gates"][:split]
    result = np.eye(16, dtype=complex)
    for gate in gates:
        if gate["gate"] == "cx":
            matrix = cnot(gate["control"], gate["target"], 4)
        else:
            matrix = embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4)
        result = matrix @ result
    return gates, torch.tensor(result, dtype=DTYPE)


def local(angles):
    mats = []
    for q in range(4):
        a, b, c = angles[q]
        za = torch.cos(a / 2) * I2 - 1j * torch.sin(a / 2) * Z
        yb = torch.cos(b / 2) * I2 - 1j * torch.sin(b / 2) * Y
        zc = torch.cos(c / 2) * I2 - 1j * torch.sin(c / 2) * Z
        mats.append(zc @ yb @ za)
    out = mats[0]
    for mat in mats[1:]:
        out = torch.kron(out, mat)
    return out


def circuit(angles, edges, prefix):
    u = prefix
    for i, edge in enumerate(edges):
        u = CX[edge] @ local(angles[i]) @ u
    return local(angles[-1]) @ u


def loss(angles, edges, prefix):
    u = circuit(angles, edges, prefix)
    m = u @ V @ u.conj().T
    off = m - torch.diag(torch.diag(m))
    return torch.sum(abs(off) ** 2) / 16


def candidate(gates, angles, edges):
    out = list(gates)
    for i in range(5):
        for q in range(4):
            for kind, angle in zip(("rz", "ry", "rz"), angles[i, q]):
                out.append({"gate": kind, "qubit": q, "theta": float(angle)})
        if i < 4:
            out.append({"gate": "cx", "control": edges[i][0],
                        "target": edges[i][1]})
    return {"n": 4, "gates": out}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", type=int, choices=range(4), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--lr", type=float, default=0.04)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    prefix_gates, prefix = prefix_data()
    edges = TOPOLOGIES[args.topology]
    angles = torch.nn.Parameter(0.3 * torch.randn(5, 4, 3, dtype=torch.float64))
    optimizer = torch.optim.Adam([angles], lr=args.lr)
    best = float("inf")
    best_angles = None
    for step in range(args.steps):
        optimizer.zero_grad()
        objective = loss(angles, edges, prefix)
        value = float(objective.detach())
        if value < best:
            best = value
            best_angles = angles.detach().clone().numpy()
        objective.backward()
        optimizer.step()
        if step % 200 == 0:
            print(json.dumps({"step": step, "loss": value, "best": best}),
                  flush=True)
    result = candidate(prefix_gates, best_angles, edges)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"topology": edges, "best_frobenius_squared_over_16": best,
                      "output": str(args.output)}))


CX = {(c, t): torch.tensor(cnot(c, t, 4), dtype=DTYPE)
      for c in range(4) for t in range(4) if c != t}


if __name__ == "__main__":
    main()
