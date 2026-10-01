"""Fit a rank-compatible five-CNOT tail after a nonlocal output gauge.

The first nine CNOTs of fivemask topology_2 are frozen.  The original six-CX
suffix is gauged by one exact two-level output rotation, while a distinct
five-CX suffix topology gets arbitrary local SU2 correction layers.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import embedded_one_qubit, evaluate, one_qubit_matrix
from delete14_nonlocal_gauge_screen import two_level, unitary
from delete14_topology_anneal import CDTYPE, RDTYPE, CX, emit, local_layer


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)


def modified_suffix(source):
    suffix = [dict(gate) for gate in source["gates"][27:]]
    cx_positions = [i for i, gate in enumerate(suffix) if gate["gate"] == "cx"]
    if len(cx_positions) != 6:
        raise ValueError("expected six-CNOT suffix")
    original_edges = [(suffix[i]["control"], suffix[i]["target"])
                      for i in cx_positions]
    if original_edges != [(2, 0), (2, 1), (2, 3), (2, 3), (0, 1), (0, 1)]:
        raise ValueError(f"unexpected suffix edges {original_edges}")
    suffix[cx_positions[0]] = {"gate": "cx", "control": 0, "target": 3}
    suffix[cx_positions[1]] = {"gate": "cx", "control": 1, "target": 3}
    del suffix[cx_positions[-1]]
    chunks, edges = [[]], []
    for gate in suffix:
        if gate["gate"] == "cx":
            edges.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert edges == [(0, 3), (1, 3), (2, 3), (2, 3), (0, 1)]
    matrices = []
    for chunk in chunks:
        matrix = np.eye(16, dtype=complex)
        for gate in chunk:
            matrix = embedded_one_qubit(one_qubit_matrix(gate),
                                         gate["qubit"], 4) @ matrix
        matrices.append(torch.as_tensor(matrix, dtype=CDTYPE))
    return chunks, edges, matrices


def suffix_matrix(angles, edges, matrices):
    out = torch.eye(16, dtype=CDTYPE)
    for slot, matrix in enumerate(matrices):
        out = local_layer(angles[slot]) @ matrix @ out
        if slot < len(edges):
            out = CX[edges[slot]] @ out
    return out


def target_loss(angles, edges, matrices, target):
    out = suffix_matrix(angles, edges, matrices)
    overlap = torch.trace(target.conj().T @ out)
    return 1 - overlap.abs() ** 2 / 256


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--adam-steps", type=int, default=300)
    p.add_argument("--bfgs-maxiter", type=int, default=150)
    p.add_argument("--lr", type=float, default=0.04)
    p.add_argument("--sigma", type=float, default=0.02)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    args = p.parse_args()
    source = json.loads((ROOT / "fivemask15_topology_2.json").read_text())
    prefix = source["gates"][:27]
    original_suffix = source["gates"][27:]
    gauge = two_level(4, 7, math.pi / 2, 0.0)
    target = torch.as_tensor(gauge @ unitary(original_suffix), dtype=CDTYPE)
    chunks, edges, matrices = modified_suffix(source)
    rng = np.random.default_rng(args.seed)
    initial_angles = rng.normal(0, args.sigma, (6, 4, 3))
    angles = torch.nn.Parameter(torch.as_tensor(initial_angles, dtype=RDTYPE))
    adam = torch.optim.Adam([angles], lr=args.lr)
    history = []
    best_loss = float("inf")
    best_angles = initial_angles.copy()
    for step in range(args.adam_steps + 1):
        adam.zero_grad()
        objective = target_loss(angles, edges, matrices, target)
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
            raise FloatingPointError("nonfinite suffix gradient")
        adam.step()

    def scipy_objective(flat):
        x = torch.tensor(flat.reshape(6, 4, 3), dtype=RDTYPE,
                         requires_grad=True)
        objective = target_loss(x, edges, matrices, target)
        gradient = torch.autograd.grad(objective, x)[0]
        return float(objective.detach()), gradient.detach().numpy().ravel().copy()

    refined = minimize(scipy_objective, best_angles.ravel(), method="BFGS",
                       jac=True, options={"maxiter": args.bfgs_maxiter,
                                          "gtol": 1e-10})
    if refined.fun < best_loss:
        best_loss, best_angles = float(refined.fun), refined.x.reshape(6, 4, 3)
    emitted_suffix = emit(chunks, edges, best_angles)["gates"]
    candidate = {"n": 4, "gates": prefix + emitted_suffix}
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    emitted_u = unitary(emitted_suffix)
    emitted_target_loss = float(1 - abs(np.trace(
        (gauge @ unitary(original_suffix)).conj().T @ emitted_u)) ** 2 / 256)
    result = {"gauge_row_pair": [4, 7], "gauge_angle_over_pi": 0.5,
              "gauge_phase_over_pi": 0, "prefix_cnot": 9,
              "suffix_edges": edges, "seed": args.seed,
              "adam_steps_cap": args.adam_steps,
              "bfgs_maxiter": args.bfgs_maxiter,
              "bfgs_nit": int(refined.nit), "bfgs_nfev": int(refined.nfev),
              "best_target_loss": best_loss,
              "emitted_gate_list_target_loss": emitted_target_loss,
              "off_diagonal_error": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"],
              "history": history}
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "history"},
                     indent=2))


if __name__ == "__main__":
    main()
