"""Dynamic positive-row weighting on the fresh chain13 numerical lead.

Three local starts, three weight-update stages each, then ordinary-cycle
polish. There are no inherited fixed gates or eigenvalue-label constraints.
All 168 axis-angle coordinates vary; every stage uses fixed positive
weights in [1,10], computed from the current output-row residuals.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, kron_gate, u3_to_axis
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_fresh_chain_refine_trial1.json"
PRIOR = ROOT / "search13_fresh_chain_refine_result.json"
OUT = ROOT / "search13_fresh_chain_dynamic_rows_result.json"
torch.set_num_threads(1)


def decode(candidate):
    chunks = [[]]
    for g in candidate["gates"]:
        if g["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(g)
    assert len(chunks) == 14
    assert all(len(chunk) == 4 and
               [g["gate"] for g in chunk] == ["u3"] * 4 and
               [g["qubit"] for g in chunk] == list(range(4))
               for chunk in chunks)
    return np.array([[u3_to_axis(g) for g in chunk] for chunk in chunks])


def unitary(angles, cxm):
    u = torch.eye(16, dtype=torch.complex128)
    for slot in range(14):
        for q in range(4):
            u = kron_gate(axis_rotation(*angles[slot, q]), q) @ u
        if slot < 13:
            u = cxm[slot] @ u
    return u


def residual_rows(angles, cxm, v):
    with torch.no_grad():
        u = unitary(torch.as_tensor(angles, dtype=torch.float64), cxm)
        d = u @ v @ u.conj().T
        off = d - torch.diag(torch.diagonal(d))
        return (off.abs() ** 2).sum(dim=1).numpy()


def weighted_objective(flat, cxm, v, weights):
    angles = torch.tensor(np.asarray(flat).reshape(14, 4, 3),
                          dtype=torch.float64, requires_grad=True)
    u = unitary(angles, cxm)
    d = u @ v @ u.conj().T
    off = d - torch.diag(torch.diagonal(d))
    loss = ((off.abs() ** 2).sum(dim=1) * weights).sum().real / weights.sum()
    gradient = torch.autograd.grad(loss, angles)[0]
    return float(loss.detach()), gradient.detach().numpy().ravel().copy()


def weighted_fit(angles, cxm, v, weights, maxiter):
    w = torch.as_tensor(weights, dtype=torch.float64)
    initial = weighted_objective(angles.ravel(), cxm, v, w)[0]
    fitted = minimize(lambda x: weighted_objective(x, cxm, v, w),
                      angles.ravel(), jac=True, method="L-BFGS-B",
                      options={"maxiter": maxiter, "ftol": 1e-15,
                               "gtol": 1e-11, "maxls": 30})
    return fitted.x.reshape(14, 4, 3), {
        "initial_weighted_loss": initial, "weighted_loss": float(fitted.fun),
        "nit": int(fitted.nit), "nfev": int(fitted.nfev),
        "success": bool(fitted.success), "message": str(fitted.message)}


def save_checked(name, schedule, angles):
    candidate = serialize([[] for _ in range(14)], schedule, angles)
    path = ROOT / name
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=34400)
    parser.add_argument("--stages", type=int, default=3)
    parser.add_argument("--weighted-maxiter", type=int, default=120)
    parser.add_argument("--direct-maxiter", type=int, default=250)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    source_angles = decode(source)
    schedule = tuple((g["control"], g["target"])
                     for g in source["gates"] if g["gate"] == "cx")
    assert len(schedule) == 13
    cxm = matrix_for_schedule(schedule)
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip, direct_grad = objective(source_angles.ravel(), fixed, cxm, v, True)
    expected = next(row["loss"] for row in json.loads(PRIOR.read_text())["rows"]
                    if row["trial"] == 1)
    assert abs(roundtrip - expected) < 1e-10
    ones_loss, ones_grad = weighted_objective(
        source_angles.ravel(), cxm, v, torch.ones(16, dtype=torch.float64))
    assert abs(ones_loss - roundtrip) < 1e-12
    assert np.max(np.abs(ones_grad - direct_grad)) < 1e-12
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "roundtrip_error": abs(roundtrip - expected),
                      "ones_loss_error": abs(ones_loss - roundtrip),
                      "ones_gradient_max_error": float(np.max(np.abs(ones_grad - direct_grad)))}), flush=True)
    rows = []
    for trial, sigma in enumerate((0.0, 0.05, 0.15)):
        seed = args.seed + trial
        angles = source_angles + np.random.default_rng(seed).normal(0, sigma, source_angles.shape)
        initial_direct = objective(angles.ravel(), fixed, cxm, v, False)
        stages = []
        for stage in range(args.stages):
            residual = residual_rows(angles, cxm, v)
            peak = float(max(residual))
            weights = 1 + 9 * residual / max(peak, 1e-30)
            assert min(weights) >= 1 and max(weights) <= 10 + 1e-12
            angles, record = weighted_fit(angles, cxm, v, weights,
                                          args.weighted_maxiter)
            ordinary_loss = objective(angles.ravel(), fixed, cxm, v, False)
            check = save_checked(
                f"search13_fresh_chain_dynamic_rows_trial{trial}_stage{stage}.json",
                schedule, angles)
            stages.append({"stage": stage, "row_residuals_before": residual.tolist(),
                           "weights": weights.tolist(), **record,
                           "ordinary_loss_after": ordinary_loss, **check})
            print(json.dumps({"trial": trial, "stage": stage,
                              "weighted_loss": record["weighted_loss"],
                              "ordinary_loss": ordinary_loss,
                              "off_diagonal_error": check["off_diagonal_error"],
                              "valid": check["valid_diagonalizer"]}), flush=True)
        polished, direct = fit(schedule, angles, fixed, v, args.direct_maxiter)
        check = save_checked(f"search13_fresh_chain_dynamic_rows_trial{trial}_direct.json",
                             schedule, polished)
        rows.append({"trial": trial, "sigma": sigma, "seed": seed,
                     "initial_direct_loss": initial_direct, "stages": stages,
                     "direct": {**direct, **check}})
        print(json.dumps({"trial": trial, "polished_loss": direct["loss"],
                          "off_diagonal_error": check["off_diagonal_error"],
                          "valid": check["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["direct"]["off_diagonal_error"])
    result = {"source": SOURCE.name,
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "schedule": schedule, "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "ones_loss_error": abs(ones_loss - roundtrip),
              "ones_gradient_max_error": float(np.max(np.abs(ones_grad - direct_grad))),
              "seed": args.seed, "stages_per_trial": args.stages,
              "weighted_maxiter": args.weighted_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "weight_update": "w_i=1+9*current row residual/max current row residual; fixed during each stage; positive in [1,10]",
              "fixed_eigenvalue_labels": False,
              "fixed_local_gates": "none; pure local layers and CNOTs",
              "rows": rows, "best_candidate": best["direct"]["candidate"],
              "best_off_diagonal_error": best["direct"]["off_diagonal_error"],
              "scope": "three starts, three dynamic positive-weight stages each, one fixed topology; bounded numerical search, not a topology exclusion"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
