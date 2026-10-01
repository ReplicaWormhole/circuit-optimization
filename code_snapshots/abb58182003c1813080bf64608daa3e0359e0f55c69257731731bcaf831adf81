"""Bounded multi-edge 14-CNOT topology annealing with local SU2 warm starts.

Every proposal changes two CNOT edges of an existing 14-CNOT topology.  A
short Adam optimization adjusts arbitrary local SU2 layers, and a cooling
schedule chooses the next topology.  The objective is the label-free
off-diagonal Frobenius loss, not an exact-circuit certificate.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch

from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift
from delete14_search import axis_rotation, axis_to_u3


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
EDGES = [(c, t) for c in range(4) for t in range(4) if c != t]
CX = {edge: torch.as_tensor(cnot(*edge, 4), dtype=CDTYPE) for edge in EDGES}
V = torch.as_tensor(right_shift(4), dtype=CDTYPE)


def split_gates(gates):
    chunks = [[]]
    edges = []
    for gate in gates:
        if gate["gate"] == "cx":
            edges.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if len(edges) != 14:
        raise ValueError(f"expected fourteen CNOTs, found {len(edges)}")
    matrices = []
    for chunk in chunks:
        matrix = np.eye(16, dtype=complex)
        for gate in chunk:
            matrix = embedded_one_qubit(
                one_qubit_matrix(gate), gate["qubit"], 4) @ matrix
        matrices.append(torch.as_tensor(matrix, dtype=CDTYPE))
    return chunks, edges, matrices


def local_layer(angles):
    result = axis_rotation(*angles[0])
    for q in range(1, 4):
        result = torch.kron(result, axis_rotation(*angles[q]))
    return result


def loss(angles, edges, matrices):
    unitary = torch.eye(16, dtype=CDTYPE)
    for slot, matrix in enumerate(matrices):
        unitary = local_layer(angles[slot]) @ matrix @ unitary
        if slot < len(edges):
            unitary = CX[edges[slot]] @ unitary
    transformed = unitary @ V @ unitary.conj().T
    offdiag = transformed - torch.diag(torch.diagonal(transformed))
    return (offdiag.abs() ** 2).sum().real / 16


def emit(chunks, edges, angles):
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q,
                          **axis_to_u3(*angles[slot, q])})
        if slot < len(edges):
            gates.append({"gate": "cx", "control": edges[slot][0],
                          "target": edges[slot][1]})
    return {"n": 4, "gates": gates}


def propose(current, rng, seen):
    for _ in range(1000):
        indices = rng.choice(len(current), size=2, replace=False)
        trial = list(current)
        for index in indices:
            alternatives = [edge for edge in EDGES if edge != trial[index]
                            and edge != trial[index][::-1]]
            trial[index] = alternatives[int(rng.integers(len(alternatives)))]
        key = tuple(trial)
        if key not in seen and sorted(map(sorted, trial)) != sorted(map(sorted, current)):
            seen.add(key)
            return trial, [int(i) for i in indices]
    raise RuntimeError("could not generate a new two-edge topology")


def optimize_short(initial, edges, matrices, steps, lr):
    angles = torch.nn.Parameter(torch.as_tensor(initial.copy(), dtype=RDTYPE))
    optimizer = torch.optim.Adam([angles], lr=lr)
    best_loss = float("inf")
    best_angles = initial.copy()
    first_loss = None
    for step in range(steps + 1):
        optimizer.zero_grad()
        objective = loss(angles, edges, matrices)
        value = float(objective.detach())
        if first_loss is None:
            first_loss = value
        if value < best_loss:
            best_loss = value
            best_angles = angles.detach().numpy().copy()
        if step == steps:
            break
        objective.backward()
        if not torch.isfinite(angles.grad).all():
            raise FloatingPointError("nonfinite SU2 gradient")
        optimizer.step()
    return first_loss, best_loss, best_angles


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--proposals", type=int, default=25)
    p.add_argument("--steps", type=int, default=12)
    p.add_argument("--lr", type=float, default=0.03)
    p.add_argument("--sigma", type=float, default=0.03)
    p.add_argument("--temp-start", type=float, default=0.15)
    p.add_argument("--temp-end", type=float, default=0.02)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    source = json.loads(args.base.read_text())
    chunks, current_edges, matrices = split_gates(source["gates"])
    current_angles = np.zeros((15, 4, 3), dtype=float)
    current_loss = float(loss(torch.as_tensor(current_angles, dtype=RDTYPE),
                              current_edges, matrices).detach())
    best_loss = current_loss
    best_edges, best_angles = list(current_edges), current_angles.copy()
    seen = {tuple(current_edges)}
    records = []
    accepted_count = 0
    for proposal_index in range(args.proposals):
        trial_edges, changed = propose(current_edges, rng, seen)
        start = current_angles + rng.normal(0, args.sigma, current_angles.shape)
        initial, optimized, trial_angles = optimize_short(
            start, trial_edges, matrices, args.steps, args.lr)
        temperature = args.temp_start * (args.temp_end / args.temp_start) ** (
            proposal_index / max(args.proposals - 1, 1))
        accepted = optimized < current_loss or rng.random() < math.exp(
            min(0.0, (current_loss - optimized) / temperature))
        if accepted:
            current_edges, current_angles, current_loss = (
                trial_edges, trial_angles, optimized)
            accepted_count += 1
        if optimized < best_loss:
            best_edges, best_angles, best_loss = (
                list(trial_edges), trial_angles.copy(), optimized)
        records.append({"proposal": proposal_index, "changed_cx": changed,
                        "topology": trial_edges, "initial_loss": initial,
                        "optimized_loss": optimized, "temperature": temperature,
                        "accepted": accepted, "best_loss": best_loss})
    candidate = emit(chunks, best_edges, best_angles)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    summary = {"base": str(args.base), "seed": args.seed,
               "proposals": args.proposals, "steps_per_proposal": args.steps,
               "temperature_start": args.temp_start,
               "temperature_end": args.temp_end,
               "accepted": accepted_count, "best_loss": best_loss,
               "best_topology": best_edges,
               "candidate_off_diagonal_error": check["off_diagonal_error"],
               "candidate_valid": check["valid_diagonalizer"],
               "records": records}
    args.result.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
