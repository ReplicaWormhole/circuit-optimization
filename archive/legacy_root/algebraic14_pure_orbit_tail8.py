"""Numerically synthesize pure-orbit V4 eigenbasis with 6+8 CNOT topology.

Freeze exact15's six-CNOT diagonal prefix. For each of nine deletions from
its nine-CNOT tail, optimize every local ZYZ layer to the pure-orbit DFT
basis of run142, modulo independent output-row phases. Screen 100 Adam
steps for all nine, refine three best for 400 further steps. This is a
bounded local search, not an exact proof or impossibility result.
"""

import json
from pathlib import Path

import numpy as np
import torch
from qiskit.quantum_info import Operator

from algebraic14_pure_orbit_householder import matched_orbit_basis
from check_circuit import evaluate
from exact_check import exact_eigenvalue_labels
from qiskit_crosscheck import to_qiskit
from topology14_eightcx_opt import emit, layers_from_suffix
from topology14_fourcx_opt import DTYPE, circuit


ROOT = Path(__file__).parent
torch.set_num_threads(1)
SCREEN_STEPS = 100
REFINE_STEPS = 400
REFINE_TOP = 3


def delete_cx(gates, deletion):
    index = -1
    result = []
    for gate in gates:
        if gate["gate"] == "cx":
            index += 1
            if index == deletion:
                continue
        result.append(gate)
    assert index == 8
    return result


def optimize(initial, edges, prefix, target, steps, seed):
    torch.manual_seed(seed)
    angles = torch.nn.Parameter(torch.tensor(initial, dtype=torch.float64))
    optimizer = torch.optim.Adam([angles], lr=0.03)
    best = float("inf")
    best_angles = None
    for _ in range(steps):
        optimizer.zero_grad()
        u = circuit(angles, edges, prefix)
        overlap = torch.diag(u @ target.conj().T)
        objective = 1 - torch.mean(torch.abs(overlap) ** 2)
        value = float(objective.detach())
        if value < best:
            best = value
            best_angles = angles.detach().clone().numpy()
        objective.backward()
        optimizer.step()
    return best, best_angles


def main():
    source = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    labels = exact_eigenvalue_labels(source)["output_labels"]
    full_source = Operator(to_qiskit(source)).data
    target_np, rows, assignments = matched_orbit_basis(
        full_source, labels)
    prefix_gates = source["gates"][:13]
    prefix_np = Operator(to_qiskit(
        {"n": 4, "gates": prefix_gates})).data
    prefix = torch.tensor(prefix_np, dtype=DTYPE)
    target = torch.tensor(target_np, dtype=DTYPE)
    suffix = source["gates"][13:]
    rows_out = []
    for deletion in range(9):
        modified = delete_cx(suffix, deletion)
        initial, edges = layers_from_suffix(modified)
        assert len(edges) == 8
        screen_loss, parameters = optimize(
            initial, edges, prefix, target, SCREEN_STEPS, deletion)
        rows_out.append({"deletion": deletion, "edges": edges,
                         "screen_loss": screen_loss,
                         "best_loss": screen_loss,
                         "parameters": parameters})
    rows_out.sort(key=lambda row: row["screen_loss"])
    for row in rows_out[:REFINE_TOP]:
        score, parameters = optimize(
            row["parameters"], row["edges"], prefix, target,
            REFINE_STEPS, row["deletion"] + 100)
        if score < row["best_loss"]:
            row["best_loss"] = score
            row["parameters"] = parameters
    rows_out.sort(key=lambda row: row["best_loss"])
    best = rows_out[0]
    candidate = emit(prefix_gates, best["parameters"], best["edges"])
    candidate_path = ROOT / "algebraic14_pure_orbit_tail8_best.json"
    candidate_path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    result = {
        "target": "pure orbit DFT matched to exact15 output labels",
        "screen_steps_per_topology": SCREEN_STEPS,
        "refine_steps": REFINE_STEPS,
        "refine_top": REFINE_TOP,
        "trials": [{key: row[key] for key in
                    ("deletion", "edges", "screen_loss", "best_loss")}
                   for row in rows_out],
        "best_deletion": best["deletion"],
        "best_target_row_phase_loss": best["best_loss"],
        "best_candidate_offdiag": check["off_diagonal_error"],
        "best_candidate_valid": check["valid_diagonalizer"],
        "best_candidate_cx": check["cnot_count"],
        "candidate": str(candidate_path),
    }
    (ROOT / "algebraic14_pure_orbit_tail8_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
