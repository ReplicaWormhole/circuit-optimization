"""Bounded seven-CNOT synthesis from an exact-rank-compatible complex gauge.

Fits one chronological schedule first to the fixed phased tail, then to the
right-cycle diagonalization equation with output eigenspace gauge free.
"""

import argparse
import json

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate, u3_to_axis
from search13_eigenrow_gauge_seven_fit import fit
from search13_joint_gauge_fit import BASE, ROOT


torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
PAIRS = ((4, 7), (8, 11))
EDGES = ((0, 2), (0, 1), (0, 3), (0, 1), (1, 2), (0, 3), (0, 1))


def target_and_prefix():
    gates = json.loads(BASE.read_text())["gates"]
    prefix, tail = gates[:13], gates[13:]
    target = np.eye(16, dtype=complex)
    for gate in tail:
        target = fixed_matrix(gate) @ target
    for first, second in PAIRS:
        a, b = target[first].copy(), target[second].copy()
        # angle pi/4 and phase -pi/2
        ep = -1j
        target[first] = (a-np.conj(ep)*b)/np.sqrt(2)
        target[second] = (ep*a+b)/np.sqrt(2)
    return prefix, target


def direct_fit(prefix, start_tail, maxiter):
    local = [g for g in start_tail if g["gate"] == "u3"]
    x0 = np.array([u3_to_axis(g) for g in local]).reshape(8, 4, 3)
    pu = np.eye(16, dtype=complex)
    for gate in prefix:
        pu = fixed_matrix(gate) @ pu
    shifted = torch.as_tensor(pu @ right_shift(4) @ pu.conj().T, dtype=CDTYPE)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                         "target": b}), dtype=CDTYPE) for a, b in EDGES]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(8):
            for qubit in range(4):
                u = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ u
            if slot < 7:
                u = cxs[slot] @ u
        transformed = u @ shifted @ u.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs()**2).sum().real/16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(x0.ravel())[0]
    opt = minimize(objective, x0.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-16,
                            "gtol": 1e-12, "maxls": 40})
    angles = opt.x.reshape(8, 4, 3)
    gates = []
    for slot in range(8):
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 7:
            a, b = EDGES[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return gates, {"initial_loss": initial, "loss": float(opt.fun),
                   "nit": int(opt.nit), "nfev": int(opt.nfev),
                   "optimizer_success": bool(opt.success),
                   "message": str(opt.message)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=25020)
    parser.add_argument("--process-maxiter", type=int, default=350)
    parser.add_argument("--direct-maxiter", type=int, default=500)
    args = parser.parse_args()
    prefix, target = target_and_prefix()
    tail, process_row = fit(EDGES, args.seed, args.process_maxiter, target)
    process_candidate = {"n": 4, "gates": prefix+tail}
    process_path = ROOT / "search13_joint_gauge_phased_process_candidate.json"
    process_path.write_text(json.dumps(process_candidate, indent=2)+"\n")
    process_check = evaluate(process_candidate)
    direct_tail, direct_row = direct_fit(prefix, tail, args.direct_maxiter)
    direct_candidate = {"n": 4, "gates": prefix+direct_tail}
    direct_path = ROOT / "search13_joint_gauge_phased_direct_candidate.json"
    direct_path.write_text(json.dumps(direct_candidate, indent=2)+"\n")
    direct_check = evaluate(direct_candidate)
    result = {"gauge_pairs": PAIRS, "angle": "pi/4",
              "phases": ["-pi/2", "-pi/2"], "edges": EDGES,
              "seed": args.seed, "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "process": {**process_row, "candidate": process_path.name,
                          "off_diagonal_error": process_check["off_diagonal_error"],
                          "valid_diagonalizer": process_check["valid_diagonalizer"]},
              "direct": {**direct_row, "candidate": direct_path.name,
                         "off_diagonal_error": direct_check["off_diagonal_error"],
                         "valid_diagonalizer": direct_check["valid_diagonalizer"]}}
    (ROOT / "search13_joint_gauge_phased_fit_result.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
