"""Two bounded seven-CNOT tails for the numerical rank-eight output gauge.

Each topology is fitted first to the gauged eight-CNOT tail by process
overlap, then continued with a label-free cycle diagonalization loss.
The gauge rank drop is numerical and the fits do not prove exact synthesis.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import torch

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate, u3_to_axis
from search13_eigenrow_gauge_seven_fit import fit as process_fit
from search13_fulltail_invariant_gauge_rank import CUTS, gauge, realign, setup

ROOT = Path(__file__).resolve().parent
GAUGE_DATA = ROOT / "search13_fulltail_invariant_gauge_rank_result.json"
OUT = ROOT / "search13_fulltail_invariant_fit_result.json"
SCHEDULES = (
    ((0, 2), (0, 3), (0, 3), (1, 3), (0, 1), (1, 2), (1, 2)),
    ((0, 2), (0, 3), (0, 3), (0, 1), (1, 2), (1, 2), (2, 3)),
)
torch.set_num_threads(1)


def decode_local_layers(tail_gates):
    chunks = [[]]
    for gate in tail_gates:
        if gate["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 8
    return np.array([[u3_to_axis(gate) for gate in chunk]
                     for chunk in chunks], dtype=float)


def direct_fit(edges, prefix, initial_angles, maxiter):
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                         "target": b}), dtype=torch.complex128)
           for a, b in edges]
    p = torch.as_tensor(prefix, dtype=torch.complex128)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=torch.float64,
                              requires_grad=True)
        u = p
        for slot in range(8):
            for qubit in range(4):
                u = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ u
            if slot < 7:
                u = cxs[slot] @ u
        diagonalized = u @ v @ u.conj().T
        off = diagonalized - torch.diag(torch.diagonal(diagonalized))
        loss = (off.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(initial_angles.ravel())[0]
    fit = minimize(objective, initial_angles.ravel(), method="L-BFGS-B",
                   jac=True, options={"maxiter": maxiter, "ftol": 1e-15,
                                      "gtol": 1e-11, "maxls": 30})
    angles = fit.x.reshape(8, 4, 3)
    gates = []
    for slot in range(8):
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 7:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return gates, {"initial_loss": initial, "loss": float(fit.fun),
                   "iterations": int(fit.nit), "success": bool(fit.success),
                   "message": str(fit.message)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process-maxiter", type=int, default=200)
    parser.add_argument("--direct-maxiter", type=int, default=300)
    parser.add_argument("--seed", type=int, default=25150)
    args = parser.parse_args()
    tail, prefix_matrix, _, sectors = setup()
    data = json.loads(GAUGE_DATA.read_text())
    best = min((row for row in data["rows"] if row["cut"] == [0, 3]),
               key=lambda row: row["loss"])
    g = gauge(torch.tensor(best["angles"], dtype=torch.float64), sectors).detach().numpy()
    target = g @ tail
    prefix_gates = json.loads((ROOT / data["source"]).read_text())["gates"][:13]
    rows = []
    for case, edges in enumerate(SCHEDULES):
        crossings = [sum((a in cut) != (b in cut) for a, b in edges)
                     for cut in CUTS]
        assert crossings[2] == 3 and min(crossings[:2]) >= 4
        tail_gates, proc = process_fit(edges, args.seed + case,
                                       args.process_maxiter, target)
        process_candidate = {"n": 4, "gates": prefix_gates + tail_gates}
        process_path = ROOT / f"search13_fulltail_invariant_fit_case{case}_process.json"
        process_path.write_text(json.dumps(process_candidate, indent=2) + "\n")
        process_check = evaluate(process_candidate)
        angles = decode_local_layers(tail_gates)
        direct_gates, direct = direct_fit(edges, prefix_matrix, angles,
                                         args.direct_maxiter)
        direct_candidate = {"n": 4, "gates": prefix_gates + direct_gates}
        direct_path = ROOT / f"search13_fulltail_invariant_fit_case{case}_direct.json"
        direct_path.write_text(json.dumps(direct_candidate, indent=2) + "\n")
        direct_check = evaluate(direct_candidate)
        rows.append({"case": case, "edges": edges, "crossing_counts": crossings,
                     "process": {**proc, "candidate": process_path.name,
                                 "off_diagonal_error": process_check["off_diagonal_error"],
                                 "valid_diagonalizer": process_check["valid_diagonalizer"]},
                     "direct": {**direct, "candidate": direct_path.name,
                                "off_diagonal_error": direct_check["off_diagonal_error"],
                                "valid_diagonalizer": direct_check["valid_diagonalizer"]}})
        print(json.dumps({"case": case, "crossings": crossings,
                          "process_loss": proc["loss"], "direct_loss": direct["loss"],
                          "direct_offdiag": direct_check["off_diagonal_error"]}),
              flush=True)
    result = {"gauge_source": GAUGE_DATA.name, "gauge_seed": best["seed"],
              "gauge_cut03_loss": best["loss"], "schedules": SCHEDULES,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter, "seed": args.seed,
              "rows": rows,
              "scope": "two graph/order fits; approximate rank-eight target; numerical only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
