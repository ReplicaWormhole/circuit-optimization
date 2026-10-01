"""Bounded identity-gauge seven-CNOT tail fits after the exact six-CNOT prefix.

Select up to six schedules passing the exact all-cut necessary screen. For each,
fit arbitrary local SU(2) layers to the fixed exact14 tail T using a random
start and a warm start obtained by deleting the nearest exact-tail CNOT.
Continue both solutions with the label-free cycle objective. Numerical fits
do not certify impossibility if they fail.
"""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, one_qubit_matrix
from delete14_search import axis_rotation, axis_to_u3, kron_gate, fixed_matrix
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
TOPOLOGY = ROOT / "search13_identity_tail_topology_result.json"
OUT = ROOT / "search13_identity_tail_fit_result.json"
torch.set_num_threads(1)


def original_tail():
    gates = json.loads(BASE.read_text())["gates"]
    prefix, tail = gates[:13], gates[13:]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    assert sum(g["gate"] == "cx" for g in tail) == 8
    pairs = tuple(tuple(sorted((g["control"], g["target"])))
                  for g in tail if g["gate"] == "cx")
    return prefix, tail, pairs


def choose(rows, count):
    groups = {}
    for row in rows:
        groups.setdefault(row["distance_to_one_exact8_deletion"], []).append(row)
    selected = []
    # The four nearest survivors are the most useful for deleting one exact
    # tail CNOT and retuning two other pair slots. Include all four before
    # broadening to more distant orders.
    for distance in sorted(groups):
        while len(selected) < count and (distance == min(groups) or
                                         (distance == min(groups) + 1 and
                                          len(selected) < count)):
            available = [row for row in groups[distance]
                         if tuple(map(tuple, row["schedule"])) not in selected]
            if not available:
                break

            def score(row):
                schedule = tuple(map(tuple, row["schedule"]))
                nearest = min((sum(a != b for a, b in zip(schedule, old))
                               for old in selected), default=7)
                return (nearest, len(set(schedule)), -row["adjacent_repeats"], schedule)

            chosen = max(available, key=score)
            selected.append(tuple(map(tuple, chosen["schedule"])))
    while len(selected) < min(count, len(rows)):
        remaining = (r for r in rows if tuple(map(tuple, r["schedule"])) not in selected)
        try:
            chosen = max(remaining, key=lambda row: (
                min(sum(a != b for a, b in zip(map(tuple, row["schedule"]), old))
                    for old in selected),
                -row["distance_to_one_exact8_deletion"],
                tuple(map(tuple, row["schedule"]))))
        except ValueError:
            break
        selected.append(tuple(map(tuple, chosen["schedule"])))
    return selected, {str(k): len(v) for k, v in sorted(groups.items())}


def directed(schedule, case):
    # Every exact-tail CNOT has control on the lower-numbered wire. Preserve
    # this convention so the deletion warm start retains its local matrices.
    return schedule


def matrix_to_axis(matrix):
    # All exact-tail one-qubit gates lie in SU(2). Remove roundoff and map
    # cos(theta/2) I - i sin(theta/2) n.sigma to its axis-angle vector.
    determinant = np.linalg.det(matrix)
    matrix = matrix / np.sqrt(determinant)
    c = np.trace(matrix).real / 2
    vector = np.array([-(matrix[0, 1] + matrix[1, 0]).imag / 2,
                       (matrix[1, 0] - matrix[0, 1]).real / 2,
                       -(matrix[0, 0] - matrix[1, 1]).imag / 2])
    s = np.linalg.norm(vector)
    if s < 1e-13:
        return np.zeros(3)
    return vector * (2 * np.arctan2(s, c) / s)


def deletion_warm(tail_gates, original_pairs, schedule):
    deletions = [original_pairs[:j] + original_pairs[j+1:]
                 for j in range(8)]
    deleted = min(range(8), key=lambda j: (
        sum(a != b for a, b in zip(schedule, deletions[j])), j))
    chunks = [[]]
    index = 0
    for gate in tail_gates:
        if gate["gate"] == "cx":
            if index != deleted:
                chunks.append([])
            index += 1
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 8 and index == 8
    angles = np.zeros((8, 4, 3))
    for slot, chunk in enumerate(chunks):
        for qubit in range(4):
            local = np.eye(2, dtype=complex)
            for gate in chunk:
                if gate["qubit"] == qubit:
                    local = one_qubit_matrix(gate) @ local
            angles[slot, qubit] = matrix_to_axis(local)
    return angles, deleted, sum(a != b for a, b in zip(schedule, deletions[deleted]))


def process_fit(edges, initial_angles, target, maxiter):
    target_t = torch.as_tensor(target, dtype=torch.complex128)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": c,
                                         "target": t}), dtype=torch.complex128)
           for c, t in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=torch.float64,
                              requires_grad=True)
        unitary = torch.eye(16, dtype=torch.complex128)
        for slot in range(8):
            for qubit in range(4):
                unitary = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ unitary
            if slot < 7:
                unitary = cxs[slot] @ unitary
        overlap = torch.trace(target_t.conj().T @ unitary)
        loss = 1 - overlap.abs().square().real / 256
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(initial_angles.ravel())[0]
    optimum = minimize(objective, initial_angles.ravel(), jac=True,
                       method="L-BFGS-B", options={"maxiter": maxiter,
                                                    "ftol": 1e-15,
                                                    "gtol": 1e-11,
                                                    "maxls": 30})
    angles = optimum.x.reshape(8, 4, 3)
    gates = []
    for slot in range(8):
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 7:
            c, t = edges[slot]
            gates.append({"gate": "cx", "control": c, "target": t})
    return gates, {"initial_loss": initial, "loss": float(optimum.fun),
                   "iterations": int(optimum.nit),
                   "success": bool(optimum.success),
                   "message": str(optimum.message)}


def checked(name, prefix, tail):
    path = ROOT / name
    candidate = {"n": 4, "gates": prefix + tail}
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=31300)
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--process-maxiter", type=int, default=180)
    parser.add_argument("--direct-maxiter", type=int, default=240)
    args = parser.parse_args()
    topology = json.loads(TOPOLOGY.read_text())
    assert topology["mixed_boundary_masks_tested"] == 127
    schedules, distance_counts = choose(topology["surviving_schedules"], args.count)
    assert schedules, "no identity-gauge rank-compatible seven-CNOT pair schedule"
    prefix, tail_gates, original_pairs = original_tail()
    tail, prefix_matrix, _, _ = setup()
    rows = []
    for case, schedule in enumerate(schedules):
        edges = directed(schedule, case)
        warm, deleted, distance = deletion_warm(tail_gates, original_pairs, schedule)
        assert distance == min(sum(a != b for a, b in zip(schedule,
                          original_pairs[:j] + original_pairs[j+1:]))
                               for j in range(8))
        rng = np.random.default_rng(args.seed + 137*case)
        starts = (("warm", warm), ("random", rng.normal(0, 0.4, (8, 4, 3))))
        trials = []
        for kind, angles in starts:
            fitted_tail, process = process_fit(edges, angles, tail, args.process_maxiter)
            process.update(checked(f"search13_identity_tail_fit_case{case}_{kind}_process.json",
                                   prefix, fitted_tail))
            direct_tail, direct = direct_fit(edges, prefix_matrix,
                                             decode_local_layers(fitted_tail),
                                             args.direct_maxiter)
            direct.update(checked(f"search13_identity_tail_fit_case{case}_{kind}_direct.json",
                                  prefix, direct_tail))
            trials.append({"start": kind, "process": process, "direct": direct})
            print(json.dumps({"case": case, "start": kind,
                              "process_loss": process["loss"],
                              "direct_loss": direct["loss"],
                              "checker_offdiagonal": direct["off_diagonal_error"],
                              "valid": direct["valid_diagonalizer"]}), flush=True)
        rows.append({"case": case, "unoriented_schedule": schedule,
                     "directed_edges": edges, "nearest_deleted_cnot": deleted,
                     "distance_to_deletion": distance,
                     "trials": trials})
    all_trials = [(r["case"], t["start"], t["direct"])
                  for r in rows for t in r["trials"]]
    best = min(all_trials, key=lambda row: row[2]["off_diagonal_error"])
    result = {"source": BASE.name, "topology_source": TOPOLOGY.name,
              "output_gauge": "identity", "seed": args.seed,
              "count": len(rows), "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "survivor_distance_counts": distance_counts,
              "selected_schedules": schedules, "rows": rows,
              "best_checker": {"case": best[0], "start": best[1], **best[2]},
              "scope": "bounded numerical fits of exact-cut-compatible pair schedules; failures are not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
