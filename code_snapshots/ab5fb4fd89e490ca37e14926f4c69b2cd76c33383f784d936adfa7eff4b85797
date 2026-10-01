"""Fit rank-compatible 14-CNOT topologies to a distinct V4 eigenbasis.

The target is the parity-conditioned fermionic Fourier matrix.  Rows are
permuted within eigenvalue sectors to maximize overlap with the 15-CNOT seed,
and arbitrary output row phases are free.  This is a numerical fit only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import linear_sum_assignment, minimize

from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift
from delete14_topology_anneal import CDTYPE, RDTYPE, CX, emit, local_layer, split_gates


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)


def unitary(gates):
    out = np.eye(16, dtype=complex)
    for gate in gates:
        matrix = (cnot(gate["control"], gate["target"], 4)
                  if gate["gate"] == "cx" else
                  embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        out = matrix @ out
    return out


def labels(u):
    diagonal = np.diag(u @ right_shift(4) @ u.conj().T)
    roots = np.array([1, 1j, -1, -1j])
    distance = abs(diagonal[:, None] - roots[None, :])
    if np.max(np.min(distance, axis=1)) > 1e-8:
        raise ValueError("alignment source or target is not a V4 eigenbasis")
    return np.argmin(distance, axis=1)


def align_target(target, source):
    source_labels, target_labels = labels(source), labels(target)
    overlap = abs(source @ target.conj().T) ** 2
    permutation = np.empty(16, dtype=int)
    for label in range(4):
        rows = np.flatnonzero(source_labels == label)
        cols = np.flatnonzero(target_labels == label)
        if len(rows) != len(cols):
            raise ValueError("eigenvalue multiplicities differ")
        i, j = linear_sum_assignment(-overlap[np.ix_(rows, cols)])
        permutation[rows[i]] = cols[j]
    return target[permutation], permutation.tolist(), float(
        overlap[np.arange(16), permutation].sum() / 16)


def circuit_matrix(angles, edges, matrices):
    out = torch.eye(16, dtype=CDTYPE)
    for slot, matrix in enumerate(matrices):
        out = local_layer(angles[slot]) @ matrix @ out
        if slot < len(edges):
            out = CX[edges[slot]] @ out
    return out


def target_loss(angles, edges, matrices, target):
    out = circuit_matrix(angles, edges, matrices)
    row_overlap = torch.diagonal(out @ target.conj().T)
    return 1 - (row_overlap.abs() ** 2).sum().real / 16


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base14", type=Path, required=True)
    p.add_argument("--seed15", type=Path, required=True)
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--topology-result", type=Path)
    p.add_argument("--proposal", type=int)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--adam-steps", type=int, default=200)
    p.add_argument("--bfgs-maxiter", type=int, default=100)
    p.add_argument("--lr", type=float, default=0.04)
    p.add_argument("--sigma", type=float, default=0.02)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    args = p.parse_args()
    base = json.loads(args.base14.read_text())
    chunks, edges, matrices = split_gates(base["gates"])
    if args.topology_result:
        if args.proposal is None:
            raise ValueError("--proposal required with --topology-result")
        batch = json.loads(args.topology_result.read_text())
        record = next(row for row in batch["records"]
                      if row["proposal"] == args.proposal)
        edges = [tuple(edge) for edge in record["topology"]]
    elif args.proposal is not None:
        raise ValueError("--proposal requires --topology-result")
    target = np.load(args.target)
    seed_u = unitary(json.loads(args.seed15.read_text())["gates"])
    aligned, permutation, seed_overlap = align_target(target, seed_u)
    target_t = torch.as_tensor(aligned, dtype=CDTYPE)
    rng = np.random.default_rng(args.seed)
    start = rng.normal(0, args.sigma, (15, 4, 3))
    angles = torch.nn.Parameter(torch.as_tensor(start, dtype=RDTYPE))
    adam = torch.optim.Adam([angles], lr=args.lr)
    history = []
    best_loss = float("inf")
    best_angles = start.copy()
    for step in range(args.adam_steps + 1):
        adam.zero_grad()
        objective = target_loss(angles, edges, matrices, target_t)
        value = float(objective.detach())
        if value < best_loss:
            best_loss, best_angles = value, angles.detach().numpy().copy()
        if step % 25 == 0:
            history.append({"adam_step": step, "loss": value,
                            "best_loss": best_loss})
        if step == args.adam_steps or value < 1e-12:
            break
        objective.backward()
        if not torch.isfinite(angles.grad).all():
            raise FloatingPointError("nonfinite target-fit gradient")
        adam.step()

    def scipy_objective(flat):
        x = torch.tensor(flat.reshape(15, 4, 3), dtype=RDTYPE,
                         requires_grad=True)
        objective = target_loss(x, edges, matrices, target_t)
        gradient = torch.autograd.grad(objective, x)[0]
        return float(objective.detach()), gradient.detach().numpy().ravel().copy()

    refined = minimize(scipy_objective, best_angles.ravel(),
                       method="BFGS", jac=True,
                       options={"maxiter": args.bfgs_maxiter, "gtol": 1e-10})
    if refined.fun < best_loss:
        best_loss, best_angles = float(refined.fun), refined.x.reshape(15, 4, 3)
    candidate = emit(chunks, edges, best_angles)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    emitted_target_loss = float(target_loss(
        torch.as_tensor(best_angles, dtype=RDTYPE), edges, matrices,
        target_t).detach())
    emitted_u = unitary(candidate["gates"])
    emitted_gate_list_loss = float(1 - np.sum(np.abs(np.diag(
        emitted_u @ aligned.conj().T)) ** 2) / 16)
    result = {"base14": str(args.base14), "seed15": str(args.seed15),
              "target": str(args.target), "topology": edges,
              "permutation": permutation,
              "seed_row_overlap_average": seed_overlap,
              "seed": args.seed, "adam_steps_cap": args.adam_steps,
              "bfgs_maxiter": args.bfgs_maxiter,
              "bfgs_nit": int(refined.nit), "bfgs_nfev": int(refined.nfev),
              "best_target_loss": best_loss,
              "emitted_internal_target_loss": emitted_target_loss,
              "emitted_gate_list_target_loss": emitted_gate_list_loss,
              "off_diagonal_error": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"],
              "history": history}
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "history"},
                     indent=2))


if __name__ == "__main__":
    main()
