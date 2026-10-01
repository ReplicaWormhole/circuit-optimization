"""All-CNOT deletion tournament from the numerically passing run359 solution.

Delete slots 0..12, keeping 14 local layers and an identity in the deleted
entangler position during optimization. Saved native circuits have exactly
12 CNOTs; the adjacent local layers require no extra counted resource.
"""

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3]))
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, right_shift
from delete14_search import axis_to_u3
from search13_adaptive_topology import matrix_for_schedule, objective
from search13_fresh_chain_dynamic_rows import decode


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
SOURCE = REPO / "collaboration/work/root/labelswap_svd_refine_candidate.json"

OUT = ROOT / "deletion_result.json"
PLAN = ROOT / "deletion_frozen.json"
torch.set_num_threads(1)


def fit_deleted(start, fixed, entanglers, v, maxiter):
    initial = objective(start.ravel(), fixed, entanglers, v, False)
    result = minimize(lambda x: objective(x, fixed, entanglers, v, True),
                      start.ravel().copy(), jac=True, method="L-BFGS-B",
                      options={"maxiter": maxiter, "ftol": 1e-15,
                               "gtol": 1e-11, "maxls": 30})
    return result.x.reshape(14, 4, 3), {
        "initial_loss": initial, "loss": float(result.fun),
        "nit": int(result.nit), "nfev": int(result.nfev),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message)}


def save_checked(name, source_schedule, deleted, angles):
    gates = []
    for slot in range(14):
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 13 and slot != deleted:
            c, t = source_schedule[slot]
            gates.append({"gate": "cx", "control": c, "target": t})
    candidate = {"n": 4, "gates": gates}
    path = ROOT / name
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    assert check["cnot_count"] == 12 and len(gates) == 14 * 4 + 12
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short-maxiter", type=int, default=100)
    parser.add_argument("--deep-maxiter", type=int, default=500)
    parser.add_argument("--deep-count", type=int, default=3)
    args = parser.parse_args()
    assert not OUT.exists() and not PLAN.exists()
    source = json.loads(SOURCE.read_text())
    start = decode(source)
    source_schedule = tuple((g["control"], g["target"])
                            for g in source["gates"] if g["gate"] == "cx")
    assert len(source_schedule) == 13
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    source_entanglers = matrix_for_schedule(source_schedule)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(start.ravel(), fixed, source_entanglers, v, False)
    from search13_fulltail_invariant_gauge_rank import circuit_matrix
    u = circuit_matrix(source["gates"])
    d = u @ right_shift(4) @ u.conj().T
    expected = float(np.sum(abs(d-np.diag(np.diag(d)))**2)/16)
    assert evaluate(source)["valid_diagonalizer"]
    assert abs(roundtrip - expected) < 1e-10
    proposals = [{"deleted_slot": slot,
                  "deleted_edge": source_schedule[slot],
                  "native_schedule": source_schedule[:slot] + source_schedule[slot+1:]}
                 for slot in range(13)]
    assert len({tuple(row["native_schedule"]) for row in proposals}) == 13
    plan = {"source": str(SOURCE.relative_to(REPO)),
            "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "source_schedule": source_schedule,
            "excluded_scope": "none; all thirteen deletions of the numerically passing run359 source",
            "optimizer_schedules": [{"deleted_slot":slot,"entanglers":[None if j==slot else edge for j,edge in enumerate(source_schedule)]} for slot in range(13)],
            "proposals": proposals}
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "roundtrip_error": abs(roundtrip - expected),
                      "planned_short_fits": 13, "planned_deep_fits": args.deep_count}), flush=True)
    short_rows = []
    angles_by_slot = {}
    def checkpoint(short_rows,deep_rows):
        path=ROOT / "deletion_checkpoint.json"
        path.write_text(json.dumps({"short_rows":short_rows,"deep_rows":deep_rows},indent=2)+"\n")
    for proposal in proposals:
        slot = proposal["deleted_slot"]
        entanglers = list(source_entanglers)
        entanglers[slot] = torch.eye(16, dtype=torch.complex128)
        angles, record = fit_deleted(start, fixed, entanglers, v, args.short_maxiter)
        check = save_checked(f"deletion_slot{slot}_short.json",
                             source_schedule, slot, angles)
        row = {**proposal, **record, **check}
        short_rows.append(row)
        angles_by_slot[slot] = angles
        checkpoint(short_rows,[])
        print(json.dumps({"stage": "short", "deleted_slot": slot,
                          "loss": row["loss"], "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
        if row["valid_diagonalizer"]:
            break
    deep_rows = []
    if not any(row["valid_diagonalizer"] for row in short_rows):
        selected = sorted(short_rows, key=lambda row: (row["loss"], row["deleted_slot"]))[:args.deep_count]
        for rank, parent in enumerate(selected):
            slot = parent["deleted_slot"]
            entanglers = list(source_entanglers)
            entanglers[slot] = torch.eye(16, dtype=torch.complex128)
            angles, record = fit_deleted(angles_by_slot[slot], fixed,
                                         entanglers, v, args.deep_maxiter)
            check = save_checked(f"deletion_slot{slot}_deep.json",
                                 source_schedule, slot, angles)
            row = {"rank": rank, "deleted_slot": slot,
                   "deleted_edge": parent["deleted_edge"],
                   "native_schedule": parent["native_schedule"],
                   "short_candidate": parent["candidate"], "short_loss": parent["loss"],
                   **record, **check}
            deep_rows.append(row)
            checkpoint(short_rows,deep_rows)
            print(json.dumps({"stage": "deep", "deleted_slot": slot,
                              "loss": row["loss"], "off_diagonal_error": row["off_diagonal_error"],
                              "valid": row["valid_diagonalizer"]}), flush=True)
            if row["valid_diagonalizer"]:
                break
    best = min(short_rows + deep_rows, key=lambda row: row["off_diagonal_error"])
    result = {**plan, "plan_file": PLAN.name,
              "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "initialization": "same source local angles for each short fit; deep starts use respective short fitted angles; no random seeds",
              "optimizer_representation": "14 local layers, 13 entangler positions, identity at deleted position; 168 free angles",
              "native_serialization": "deleted entangler omitted, exactly 12 CNOTs, adjacent local layers retained",
              "short_maxiter": args.short_maxiter, "deep_maxiter": args.deep_maxiter,
              "short_count": len(short_rows), "deep_count": len(deep_rows),
              "short_rows": short_rows, "deep_rows": deep_rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "all thirteen deletions, one source warm start each, best three continued by fitted loss; bounded numerical failure is not a lower bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
