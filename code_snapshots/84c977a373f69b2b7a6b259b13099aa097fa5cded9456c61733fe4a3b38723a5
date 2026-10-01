"""Direct gauge-free fit of causally diverse seven-CNOT full-tail schedules.

The exact six-CNOT prefix remains fixed.  Each schedule interleaves crossing
edges rather than preserving the 02/13 and 01/23 block order.  The loss is
only the right-cycle off-diagonal Frobenius norm; output eigenspace gauge and
eigenvalue order are free.  A prior seven-CNOT candidate supplies a topology-
transfer warm start, followed by continuation from the best prior schedule.
"""

import argparse
import json

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate, u3_to_axis
from search13_joint_gauge_fit import BASE, ROOT


torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
WARM = ROOT / "search13_joint_gauge_phased_process_candidate.json"
SCHEDULES = (
    ((0, 2), (0, 1), (2, 3), (0, 2), (1, 2), (0, 1), (0, 3)),
    ((0, 3), (1, 2), (0, 1), (2, 3), (0, 2), (1, 3), (0, 3)),
    ((1, 3), (0, 2), (0, 3), (1, 2), (2, 3), (0, 1), (1, 2)),
)
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
REQUIRED_RANKS = (4, 4, 4, 4, 16, 16, 16)


def check_graph(edges):
    cap = tuple(2**sum((a in cut) != (b in cut) for a, b in edges)
                for cut in CUTS)
    assert all(r <= c for r, c in zip(REQUIRED_RANKS, cap))
    assert len(set(tuple(sorted(edge)) for edge in edges)) >= 5
    assert all(edges[i] != edges[i+1] for i in range(6))
    return cap


def prepare_prefix():
    gates = json.loads(BASE.read_text())["gates"][:13]
    assert sum(g["gate"] == "cx" for g in gates) == 6
    pu = np.eye(16, dtype=complex)
    for gate in gates:
        pu = fixed_matrix(gate) @ pu
    shifted = pu @ right_shift(4) @ pu.conj().T
    return gates, torch.as_tensor(shifted, dtype=CDTYPE)


def warm_angles():
    gates = json.loads(WARM.read_text())["gates"][13:]
    local = [g for g in gates if g["gate"] == "u3"]
    assert len(local) == 32
    return np.array([u3_to_axis(g) for g in local]).reshape(8, 4, 3)


def fit(edges, shifted, initial, maxiter):
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                         "target": b}), dtype=CDTYPE) for a, b in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(8):
            for wire in range(4):
                u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
            if slot < 7:
                u = cxs[slot] @ u
        transformed = u @ shifted @ u.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs()**2).sum().real/16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    start_loss = objective(initial.ravel())[0]
    opt = minimize(objective, initial.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-16,
                            "gtol": 1e-12, "maxls": 40})
    return opt.x.reshape(8, 4, 3), {"initial_loss": start_loss,
           "loss": float(opt.fun), "nit": int(opt.nit), "nfev": int(opt.nfev),
           "optimizer_success": bool(opt.success), "message": str(opt.message)}


def candidate(prefix, edges, angles):
    gates = list(prefix)
    for slot in range(8):
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 7:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=25030)
    parser.add_argument("--maxiter", type=int, default=350)
    parser.add_argument("--sigma", type=float, default=0.4)
    args = parser.parse_args()
    prefix, shifted = prepare_prefix()
    rng = np.random.default_rng(args.seed)
    transferred = warm_angles()
    transferred_source = WARM.name
    rows = []
    for index, edges in enumerate(SCHEDULES):
        caps = check_graph(edges)
        local_rows = []
        for start_kind, start in (("transferred", transferred),
                                  ("random", rng.normal(0, args.sigma, (8, 4, 3)))):
            angles, record = fit(edges, shifted, start.copy(), args.maxiter)
            circuit = candidate(prefix, edges, angles)
            path = ROOT / f"search13_joint_gauge_interleaved_s{index}_{start_kind}.json"
            path.write_text(json.dumps(circuit, indent=2)+"\n")
            checked = evaluate(circuit)
            record.update(schedule_index=index, start_kind=start_kind,
                          start_source=(transferred_source if start_kind == "transferred"
                                        else f"rng_seed_{args.seed}"),
                          edges=edges, cut_capacities=caps,
                          candidate=path.name,
                          off_diagonal_error=checked["off_diagonal_error"],
                          valid_diagonalizer=checked["valid_diagonalizer"])
            local_rows.append((record, angles))
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
        best_record, transferred = min(local_rows, key=lambda pair: pair[0]["loss"])
        transferred_source = best_record["candidate"]
    result = {"base": BASE.name, "fixed_prefix_gate_count": len(prefix),
              "warm_start": WARM.name, "schedules": SCHEDULES,
              "seed": args.seed, "sigma": args.sigma,
              "maxiter_per_fit": args.maxiter, "rows": rows,
              "best_loss": min(rows, key=lambda row: row["loss"]),
              "best_off_diagonal": min(rows, key=lambda row: row["off_diagonal_error"])}
    (ROOT / "search13_joint_gauge_interleaved_direct_result.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
