"""Whole-circuit 13-CNOT fit with a changed prefix and gauge-orbit stages.

Topology proposals are simultaneous mutations of the exact 14-CNOT circuit:
delete a prefix CNOT, reroute a different retained prefix CNOT, and reroute
both CNOTs of one suffix matchgate block. All 168 local SU(2) parameters
across the 13-CNOT circuit vary in every optimization stage.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, kron_gate
from search13_adaptive_topology import objective as direct_objective
from search13_fulltail_invariant_gauge_rank import circuit_matrix, setup
from search13_multi_pair_block_rewire import prepare
from search13_procrustes_gauge_fit import optimal_gauge


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_joint_prefix_orbit_fit_result.json"
CASES = (
    # deleted old prefix CNOT, changed retained prefix CNOT, four suffix blocks
    (4, (1, 0, 3), ((0, 2), (1, 3), (0, 3), (2, 3))),
    (5, (2, 0, 1), ((0, 2), (0, 3), (0, 1), (2, 3))),
)
torch.set_num_threads(1)


def unitary(angles, fixed, cxs):
    u = torch.eye(16, dtype=torch.complex128)
    for slot in range(14):
        u = fixed[slot] @ u
        for wire in range(4):
            u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
        if slot < 13:
            u = cxs[slot] @ u
    return u


def process_objective(flat, fixed, cxs, target):
    angles = torch.tensor(np.asarray(flat).reshape(14, 4, 3),
                          dtype=torch.float64, requires_grad=True)
    u = unitary(angles, fixed, cxs)
    overlap = torch.trace(target.conj().T @ u)
    loss = 1 - overlap.abs().square().real / 256
    grad = torch.autograd.grad(loss, angles)[0]
    return float(loss.detach()), grad.detach().numpy().ravel().copy()


def fit(objective, start, maxiter):
    result = minimize(objective, np.asarray(start).ravel(), jac=True,
                      method="L-BFGS-B", options={"maxiter": maxiter,
                                                  "ftol": 1e-15,
                                                  "gtol": 1e-11,
                                                  "maxls": 30})
    return result.x.reshape(14, 4, 3), {
        "loss": float(result.fun), "iterations": int(result.nit),
        "success": bool(result.success), "message": str(result.message)}


def serialize(chunks, gates, angles):
    out = []
    for slot, chunk in enumerate(chunks):
        out.extend(chunk)
        for wire in range(4):
            out.append({"gate": "u3", "qubit": wire,
                        **axis_to_u3(*angles[slot, wire])})
        if slot < 13:
            out.append(gates[slot])
    return {"n": 4, "gates": out}


def checked(path, candidate):
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    result = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": result["cnot_count"],
            "off_diagonal_error": result["off_diagonal_error"],
            "valid_diagonalizer": result["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=32700)
    parser.add_argument("--outer", type=int, default=2)
    parser.add_argument("--process-maxiter", type=int, default=70)
    parser.add_argument("--direct-maxiter", type=int, default=180)
    args = parser.parse_args()
    source = json.loads(BASE.read_text())
    exact = circuit_matrix(source["gates"])
    _, _, labels, sectors = setup()
    rows = []
    for case, (deleted, rewire, blocks) in enumerate(CASES):
        chunks, cxs, fixed, cxm = prepare(deleted, blocks,
                                            prefix_rewire=rewire)
        for start_name in ("zero", "random"):
            rng = np.random.default_rng(args.seed + 11 * case +
                                        (start_name == "random"))
            angles = (np.zeros((14, 4, 3)) if start_name == "zero"
                      else rng.normal(0, 0.15, (14, 4, 3)))
            stages = []
            for stage in range(args.outer):
                with torch.no_grad():
                    current = unitary(torch.as_tensor(angles), fixed,
                                      cxm).numpy()
                gauge = optimal_gauge(current, exact, sectors)
                assert np.max(np.abs(gauge @ gauge.conj().T - np.eye(16))) < 1e-12
                assert all(abs(gauge[i, j]) < 1e-12 or labels[i] == labels[j]
                           for i in range(16) for j in range(16))
                target = torch.as_tensor(gauge @ exact, dtype=torch.complex128)
                angles, record = fit(
                    lambda x: process_objective(x, fixed, cxm, target),
                    angles, args.process_maxiter)
                stages.append(record)
            process_candidate = serialize(chunks, cxs, angles)
            process_path = ROOT / f"search13_joint_prefix_orbit_case{case}_{start_name}_process.json"
            process_check = checked(process_path, process_candidate)
            v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
            angles, direct = fit(
                lambda x: direct_objective(x, fixed, cxm, v, True),
                angles, args.direct_maxiter)
            direct_candidate = serialize(chunks, cxs, angles)
            direct_path = ROOT / f"search13_joint_prefix_orbit_case{case}_{start_name}_direct.json"
            direct_check = checked(direct_path, direct_candidate)
            row = {"case": case, "deleted_prefix_cnot": deleted,
                   "prefix_rewire": rewire, "suffix_blocks": blocks,
                   "start": start_name, "stages": stages,
                   "process_check": process_check,
                   "direct": {**direct, **direct_check}}
            rows.append(row)
            print(json.dumps({"case": case, "start": start_name,
                              "process_loss": stages[-1]["loss"],
                              "direct_loss": direct["loss"],
                              "off_diagonal_error": direct_check["off_diagonal_error"],
                              "valid": direct_check["valid_diagonalizer"]}),
                  flush=True)
    best = min(rows, key=lambda r: r["direct"]["off_diagonal_error"])
    result = {"source": BASE.name,
              "strategy": "simultaneous prefix deletion, distant prefix rewire, and suffix-block edge change; full local layers fitted to optimal commuting output gauge, then label-free loss",
              "cases": CASES, "seed": args.seed, "outer": args.outer,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "rows": rows,
              "best_candidate": best["direct"]["candidate"],
              "best_off_diagonal_error": best["direct"]["off_diagonal_error"],
              "scope": "four bounded nonconvex fits; no lower-bound implication"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
