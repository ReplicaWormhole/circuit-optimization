"""Gauge-free continuation of the best interleaved seven-CNOT tail fit."""

import argparse
import json

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate, u3_to_axis
from search13_joint_gauge_fit import ROOT


torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
SOURCE = ROOT / "search13_joint_gauge_candidate_1.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=500)
    args = parser.parse_args()
    source = json.loads(SOURCE.read_text())
    prefix, tail = source["gates"][:13], source["gates"][13:]
    edges = [(g["control"], g["target"]) for g in tail if g["gate"] == "cx"]
    assert len(edges) == 7 and len(tail) == 39
    local = [g for g in tail if g["gate"] == "u3"]
    assert len(local) == 32
    x0 = np.array([u3_to_axis(g) for g in local]).reshape(8, 4, 3)
    pu = np.eye(16, dtype=complex)
    for gate in prefix:
        pu = fixed_matrix(gate) @ pu
    shifted = torch.as_tensor(pu @ right_shift(4) @ pu.conj().T, dtype=CDTYPE)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                         "target": b}), dtype=CDTYPE) for a, b in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(8, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        unitary = torch.eye(16, dtype=CDTYPE)
        for slot in range(8):
            for qubit in range(4):
                unitary = kron_gate(axis_rotation(*angles[slot, qubit]), qubit) @ unitary
            if slot < 7:
                unitary = cxs[slot] @ unitary
        transformed = unitary @ shifted @ unitary.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs()**2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial = objective(x0.ravel())[0]
    fit = minimize(objective, x0.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": args.maxiter, "ftol": 1e-16,
                            "gtol": 1e-12, "maxls": 40})
    angles = fit.x.reshape(8, 4, 3)
    gates = list(prefix)
    for slot in range(8):
        for qubit in range(4):
            gates.append({"gate": "u3", "qubit": qubit,
                          **axis_to_u3(*angles[slot, qubit])})
        if slot < 7:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    candidate = {"n": 4, "gates": gates}
    path = ROOT / "search13_joint_gauge_continuation_candidate.json"
    path.write_text(json.dumps(candidate, indent=2)+"\n")
    check = evaluate(candidate)
    result = {"source": SOURCE.name, "edges": edges, "maxiter": args.maxiter,
              "initial_loss": initial, "final_loss": float(fit.fun),
              "nit": int(fit.nit), "nfev": int(fit.nfev),
              "optimizer_success": bool(fit.success), "message": str(fit.message),
              "candidate": path.name,
              "off_diagonal_error": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"]}
    (ROOT / "search13_joint_gauge_continuation_result.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
