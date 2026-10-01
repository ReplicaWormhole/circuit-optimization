"""Move a numerical 14-CNOT diagonalizer toward a rational-pi angle grid.

At a valid gate list, numerically estimate the Jacobian of U V4 - D U.
Project a grid-penalty gradient into its numerical nullspace, then correct
back to the diagonalizer manifold by Gauss-Newton.  This is a bounded gauge
search.  A small residual and a numerical Jacobian rank are not exact proof.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import embedded_one_qubit, evaluate, one_qubit_matrix
from delete14_nearmiss_refine import (ROOTS, V, emit, parse, residual, unitary)
from check_circuit import cnot


def gate_list_unitary(gates):
    out = np.eye(16, dtype=complex)
    for gate in gates:
        matrix = (cnot(gate["control"], gate["target"], 4)
                  if gate["gate"] == "cx" else
                  embedded_one_qubit(one_qubit_matrix(gate),
                                     gate["qubit"], 4))
        out = matrix @ out
    return out


def mixer_magnitudes(u, reference):
    w = u @ reference.conj().T
    return {"5_5": float(abs(w[5, 5])), "5_9": float(abs(w[5, 9])),
            "6_6": float(abs(w[6, 6])), "6_10": float(abs(w[6, 10]))}


def penalty(angles, denominator):
    return float(np.sum(np.sin(denominator * angles) ** 2))


def grid_counts(angles, denominator):
    step = math.pi / denominator
    distance = abs(angles - np.round(angles / step) * step)
    return {str(tolerance): int(np.count_nonzero(distance < tolerance))
            for tolerance in (1e-6, 1e-4, 1e-3, 1e-2)}


def jacobian(angles, prefix_u, cx_matrices, labels, difference_step):
    flat = angles.ravel()
    base = residual(unitary(prefix_u, cx_matrices, angles), labels)
    matrix = np.empty((len(base), len(flat)), dtype=float)
    for coordinate in range(len(flat)):
        plus, minus = flat.copy(), flat.copy()
        plus[coordinate] += difference_step
        minus[coordinate] -= difference_step
        rp = residual(unitary(prefix_u, cx_matrices,
                              plus.reshape(9, 4, 3)), labels)
        rm = residual(unitary(prefix_u, cx_matrices,
                              minus.reshape(9, 4, 3)), labels)
        matrix[:, coordinate] = (rp - rm) / (2 * difference_step)
    return matrix, base


def correct(angles, prefix_u, cx_matrices, labels, difference_step,
            rcond, rounds):
    trial = angles.copy()
    for _ in range(rounds):
        matrix, r = jacobian(trial, prefix_u, cx_matrices, labels,
                             difference_step)
        if np.linalg.norm(r) < 1e-12:
            break
        delta = np.linalg.lstsq(matrix, -r, rcond=rcond)[0]
        trial = (trial.ravel() + delta).reshape(9, 4, 3)
    return trial


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--grid-denominator", type=int, default=16)
    p.add_argument("--maxiter", type=int, default=20)
    p.add_argument("--step-norm", type=float, default=0.03)
    p.add_argument("--difference-step", type=float, default=1e-5)
    p.add_argument("--rcond", type=float, default=1e-9)
    args = p.parse_args()
    source = json.loads(args.input.read_text())
    prefix, prefix_u, edges, angles = parse(source)
    cx_matrices = [cnot(*edge, 4) for edge in edges]
    diagonal = np.diag(unitary(prefix_u, cx_matrices, angles)
                       @ V @ unitary(prefix_u, cx_matrices, angles).conj().T)
    labels = ROOTS[np.argmin(abs(diagonal[:, None] - ROOTS[None, :]), axis=1)]
    reference = gate_list_unitary(json.loads((
        Path(__file__).resolve().parent /
        "topology16_15_exact_candidate.json").read_text())["gates"])
    initial_mixer = mixer_magnitudes(unitary(prefix_u, cx_matrices, angles),
                                     reference)
    initial_penalty = penalty(angles, args.grid_denominator)
    initial_counts = grid_counts(angles, args.grid_denominator)
    history = []
    for iteration in range(args.maxiter):
        matrix, r = jacobian(angles, prefix_u, cx_matrices, labels,
                             args.difference_step)
        _, singular, vt = np.linalg.svd(matrix, full_matrices=False)
        rank = int(np.count_nonzero(singular > args.rcond * singular[0]))
        nullspace = vt[rank:].T
        gradient = args.grid_denominator * np.sin(
            2 * args.grid_denominator * angles.ravel())
        direction = -nullspace @ (nullspace.T @ gradient)
        norm = float(np.linalg.norm(direction))
        if norm < 1e-10:
            history.append({"iteration": iteration, "stopped": "flat-null-gradient",
                            "jacobian_rank": rank})
            break
        direction *= args.step_norm / norm
        old_penalty = penalty(angles, args.grid_denominator)
        accepted = False
        for factor in (1.0, 0.5, 0.25, 0.125):
            trial = (angles.ravel() + factor * direction).reshape(9, 4, 3)
            trial = correct(trial, prefix_u, cx_matrices, labels,
                            args.difference_step, args.rcond, 2)
            trial_residual = float(np.linalg.norm(residual(
                unitary(prefix_u, cx_matrices, trial), labels)))
            trial_penalty = penalty(trial, args.grid_denominator)
            if trial_residual < 1e-9 and trial_penalty < old_penalty:
                angles = trial
                accepted = True
                break
        history.append({"iteration": iteration, "jacobian_rank": rank,
                        "null_dimension": len(angles.ravel()) - rank,
                        "old_penalty": old_penalty,
                        "new_penalty": trial_penalty if accepted else None,
                        "accepted_factor": factor if accepted else None,
                        "residual_norm": trial_residual if accepted else None})
        if not accepted:
            break
    candidate = emit(prefix, edges, angles)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    result = {"input": str(args.input),
              "grid_denominator": args.grid_denominator,
              "maxiter": args.maxiter, "step_norm": args.step_norm,
              "difference_step": args.difference_step, "rcond": args.rcond,
              "initial_penalty": initial_penalty,
              "final_penalty": penalty(angles, args.grid_denominator),
              "initial_grid_counts": initial_counts,
              "final_grid_counts": grid_counts(angles, args.grid_denominator),
              "initial_mixer_magnitudes": initial_mixer,
              "final_mixer_magnitudes": mixer_magnitudes(
                  unitary(prefix_u, cx_matrices, angles), reference),
              "final_offdiag": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"],
              "cnot_count": check["cnot_count"], "history": history}
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "history"},
                     indent=2))


if __name__ == "__main__":
    main()
