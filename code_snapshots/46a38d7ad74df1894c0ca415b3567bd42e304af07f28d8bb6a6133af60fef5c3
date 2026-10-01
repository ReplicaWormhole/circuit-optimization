"""Optimize all local gates of an eight-CNOT Fourier tail after six-CNOT prefix.

The tail topology keeps the two exact butterflies and the final two F2
blocks but removes the intervening CZ. Initialization is the exact 15-CNOT
gate list with that CZ omitted, canonically regrouped into ZYZ layers.
Every layer then varies; the off-diagonal objective allows a new degenerate
eigenbasis and output eigenvalue order.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from qiskit.synthesis import OneQubitEulerDecomposer

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix
from topology14_fourcx_opt import DTYPE, loss


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)
DECOMPOSER = OneQubitEulerDecomposer("ZYZ")


def source_data():
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix = source["gates"][:13]
    assert source["gates"][49:52] == [
        {"gate": "h", "qubit": 1},
        {"gate": "cx", "control": 2, "target": 1},
        {"gate": "h", "qubit": 1},
    ]
    suffix = source["gates"][13:49] + source["gates"][52:]
    u = np.eye(16, dtype=complex)
    for gate in prefix:
        m = (cnot(gate["control"], gate["target"], 4) if gate["gate"] == "cx"
             else embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        u = m @ u
    return prefix, torch.tensor(u, dtype=DTYPE), suffix


def layers_from_suffix(suffix):
    wires = [np.eye(2, dtype=complex) for _ in range(4)]
    layers = []
    edges = []
    for gate in suffix:
        if gate["gate"] != "cx":
            q = gate["qubit"]
            wires[q] = one_qubit_matrix(gate) @ wires[q]
            continue
        layers.append(wires)
        edges.append((gate["control"], gate["target"]))
        wires = [np.eye(2, dtype=complex) for _ in range(4)]
    layers.append(wires)
    out = np.zeros((len(layers), 4, 3), dtype=float)
    for i, layer in enumerate(layers):
        for q, matrix in enumerate(layer):
            theta, phi, lam = DECOMPOSER.angles(matrix)
            out[i, q] = (lam, theta, phi)
    return out, edges


def emit(prefix, parameters, edges):
    gates = list(prefix)
    for i in range(len(edges) + 1):
        for q in range(4):
            for name, theta in zip(("rz", "ry", "rz"), parameters[i, q]):
                gates.append({"gate": name, "qubit": q, "theta": float(theta)})
        if i < len(edges):
            gates.append({"gate": "cx", "control": edges[i][0],
                          "target": edges[i][1]})
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--lr", type=float, default=0.03)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    prefix, prefix_u, suffix = source_data()
    initial, edges = layers_from_suffix(suffix)
    assert len(edges) == 8
    parameters = torch.nn.Parameter(torch.tensor(initial, dtype=torch.float64) +
                                    0.01 * torch.randn(initial.shape, dtype=torch.float64))
    optimizer = torch.optim.Adam([parameters], lr=args.lr)
    initial_loss = float(loss(torch.tensor(initial, dtype=torch.float64), edges,
                              prefix_u).detach())
    best = float("inf")
    best_parameters = None
    for step in range(args.steps):
        optimizer.zero_grad()
        objective = loss(parameters, edges, prefix_u)
        value = float(objective.detach())
        if value < best:
            best = value
            best_parameters = parameters.detach().clone().numpy()
        objective.backward()
        optimizer.step()
        if step % 200 == 0:
            print(json.dumps({"step": step, "loss": value, "best": best}),
                  flush=True)
    args.output.write_text(json.dumps(emit(prefix, best_parameters, edges),
                                      indent=2) + "\n")
    print(json.dumps({"initial_loss": initial_loss, "best_loss": best,
                      "topology": edges, "output": str(args.output)}))


if __name__ == "__main__":
    main()
