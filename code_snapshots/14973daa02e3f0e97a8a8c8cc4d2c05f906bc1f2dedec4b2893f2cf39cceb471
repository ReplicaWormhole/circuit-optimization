"""Adaptive 13-CNOT topology mutations from a prefix-deletion near-hit.

The exact14 local chunks after deletion of prefix CX4 stay as fixed local
matrices; arbitrary SU(2) layers between every CNOT are refitted. A beam of
near-hit 13-CNOT schedules proposes single changed CNOT pairs anywhere in
the circuit. Raw-loss screening is followed by bounded L-BFGS-B fitting.
All results are numerical search evidence, not a lower bound.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import BLOCKS, SOURCE, decode_layers


ROOT = Path(__file__).resolve().parent
PAIRS = tuple((a, b) for a in range(4) for b in range(4) if a != b)
CDTYPE = torch.complex128
RDTYPE = torch.float64
torch.set_num_threads(1)


def schedule_key(schedule):
    return tuple(tuple(edge) for edge in schedule)


def admissible(schedule):
    deg = [sum(wire in edge for edge in schedule) for wire in range(4)]
    if min(deg) < 2:
        return False
    for wire in range(4):
        cone = {wire}
        for edge in schedule:
            if cone.intersection(edge):
                cone.update(edge)
        if len(cone) != 4:
            return False
    return True


def matrix_for_schedule(schedule):
    return [torch.as_tensor(fixed_matrix({"gate": "cx", "control": c,
                                           "target": t}), dtype=CDTYPE)
            for c, t in schedule]


def objective(flat, fixed, cxm, v, gradient):
    angles = torch.tensor(np.asarray(flat).reshape(14, 4, 3), dtype=RDTYPE,
                          requires_grad=gradient)
    u = torch.eye(16, dtype=CDTYPE)
    for slot in range(14):
        u = fixed[slot] @ u
        for wire in range(4):
            u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
        if slot < 13:
            u = cxm[slot] @ u
    d = u @ v @ u.conj().T
    off = d - torch.diag(torch.diagonal(d))
    loss = (off.abs() ** 2).sum().real / 16
    if gradient:
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()
    return float(loss.detach())


def fit(schedule, angles0, fixed, v, maxiter):
    cxm = matrix_for_schedule(schedule)
    x0 = np.asarray(angles0, dtype=float).ravel().copy()
    initial = objective(x0, fixed, cxm, v, False)
    result = minimize(lambda x: objective(x, fixed, cxm, v, True), x0,
                      jac=True, method="L-BFGS-B",
                      options={"maxiter": maxiter, "ftol": 1e-15,
                               "gtol": 1e-11, "maxls": 30})
    return result.x.reshape(14, 4, 3), {"initial_loss": initial,
              "loss": float(result.fun), "nit": int(result.nit),
              "nfev": int(result.nfev), "optimizer_success": bool(result.success),
              "optimizer_message": str(result.message)}


def serialize(chunks, schedule, angles):
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 13:
            c, t = schedule[slot]
            gates.append({"gate": "cx", "control": c, "target": t})
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generations", type=int, default=4)
    parser.add_argument("--proposals", type=int, default=24)
    parser.add_argument("--fits", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=150)
    parser.add_argument("--seed", type=int, default=135000)
    parser.add_argument("--prefix", default="search13_adaptive_topology")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    chunks, cxs, fixed, _ = prepare(4, BLOCKS)
    start_schedule = [(gate["control"], gate["target"]) for gate in cxs]
    start_angles = decode_layers(json.loads(SOURCE.read_text()))
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    start_loss = objective(start_angles.ravel(), fixed,
                           matrix_for_schedule(start_schedule), v, False)
    if abs(start_loss - 0.14176255852774977) > 1e-8:
        raise AssertionError(f"source candidate roundtrip loss {start_loss}")
    beam = [(start_loss, start_schedule, start_angles)]
    seen = {schedule_key(start_schedule)}
    records = []
    for generation in range(args.generations):
        proposals = []
        for parent_id, (_, parent_schedule, parent_angles) in enumerate(beam):
            generated = 0
            attempts = 0
            while generated < args.proposals and attempts < 20 * args.proposals:
                attempts += 1
                slot = int(rng.integers(13))
                old = parent_schedule[slot]
                pair = PAIRS[int(rng.integers(len(PAIRS)))]
                if pair == old or set(pair) == set(old):
                    continue
                trial = list(parent_schedule)
                trial[slot] = pair
                key = schedule_key(trial)
                if key in seen or not admissible(trial):
                    continue
                seen.add(key)
                generated += 1
                raw = objective(parent_angles.ravel(), fixed,
                                matrix_for_schedule(trial), v, False)
                proposals.append((raw, parent_id, slot, old, pair, trial,
                                  parent_angles))
        proposals.sort(key=lambda x: x[0])
        fitted = []
        for rank, (raw, parent_id, slot, old, pair, schedule,
                   parent_angles) in enumerate(proposals[:args.fits]):
            angles, fit_record = fit(schedule, parent_angles, fixed, v,
                                     args.maxiter)
            candidate = serialize(chunks, schedule, angles)
            check = evaluate(candidate)
            path = ROOT / f"{args.prefix}_g{generation}_r{rank}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            row = {"generation": generation, "rank": rank,
                   "parent_id": parent_id, "changed_slot": slot,
                   "old_pair": old, "new_pair": pair, "raw_loss": raw,
                   **fit_record, "candidate": path.name,
                   "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid_diagonalizer": check["valid_diagonalizer"]}
            records.append(row)
            fitted.append((fit_record["loss"], schedule, angles))
            print(json.dumps(row, sort_keys=True), flush=True)
        beam = sorted(beam + fitted, key=lambda x: x[0])[:2]
        if not fitted:
            break
    result = {"source": SOURCE.name, "seed": args.seed,
              "generations": args.generations, "proposals_per_parent": args.proposals,
              "fits_per_generation": args.fits, "maxiter": args.maxiter,
              "start_loss": start_loss, "unique_schedules_seen": len(seen),
              "records": records,
              "best": min(records, key=lambda row: row["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": path.name,
                      "best": result["best"]["candidate"]}), flush=True)


if __name__ == "__main__":
    main()
