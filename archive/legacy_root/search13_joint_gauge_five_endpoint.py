"""Joint all-local fit of new five-CNOT prefix endpoints and eight-CNOT tails.

The five-edge prefixes are wholesale new linear networks, not deletions or
single-edge rewires of the exact14 prefix.  Each exposes full input parity on
one wire.  Every local SU(2) layer of the resulting 13-CNOT circuit varies
under the label-free cycle loss, including both sides of the prefix boundary.
"""

import argparse
import json

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import evaluate, one_qubit_matrix, right_shift
from delete14_search import axis_rotation, axis_to_u3, fixed_matrix, kron_gate
from search13_joint_gauge_fit import ROOT


torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64
WARM = ROOT / "search13_prefix_delete_fit_p4_s1.json"
PREFIXES = (
    ((3, 1), (1, 2), (0, 3), (2, 0), (3, 2)),
    ((1, 3), (3, 2), (0, 1), (2, 0), (1, 2)),
)
EXACT_PREFIX = ((3, 1), (1, 2), (3, 2), (0, 2), (1, 2), (3, 2))
TAILS = (
    ((0, 2), (0, 2), (1, 3), (1, 3),
     (0, 1), (0, 1), (2, 3), (2, 3)),
    ((0, 2), (1, 3), (0, 2), (1, 3),
     (0, 1), (2, 3), (0, 1), (2, 3)),
)


def endpoint(edges):
    rows = [8, 4, 2, 1]
    for control, target in edges:
        rows[target] ^= rows[control]
    return tuple(rows)


def local_to_axis(unitary):
    # Remove an arbitrary U(2) phase, then invert the axis-angle exponential.
    u = unitary / np.sqrt(np.linalg.det(unitary))
    c = (np.trace(u)/2).real
    vector = np.array([-(u[0, 1]+u[1, 0]).imag/2,
                       (u[1, 0]-u[0, 1]).real/2,
                       -(u[0, 0]-u[1, 1]).imag/2])
    sine = np.linalg.norm(vector)
    if sine < 1e-13:
        return np.zeros(3)
    return vector * (2*np.arctan2(sine, c)/sine)


def warm_angles():
    source = json.loads(WARM.read_text())["gates"]
    chunks = [[]]
    for gate in source:
        if gate["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 14
    layers = []
    for chunk in chunks:
        mats = [np.eye(2, dtype=complex) for _ in range(4)]
        for gate in chunk:
            q = gate["qubit"]
            mats[q] = one_qubit_matrix(gate) @ mats[q]
        layers.append([local_to_axis(mat) for mat in mats])
    return np.asarray(layers)


def fit(edges, initial, maxiter):
    v = torch.as_tensor(right_shift(4), dtype=CDTYPE)
    cxs = [torch.as_tensor(fixed_matrix({"gate": "cx", "control": a,
                                         "target": b}), dtype=CDTYPE) for a, b in edges]

    def objective(flat):
        angles = torch.tensor(flat.reshape(14, 4, 3), dtype=RDTYPE,
                              requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(14):
            for wire in range(4):
                u = kron_gate(axis_rotation(*angles[slot, wire]), wire) @ u
            if slot < 13:
                u = cxs[slot] @ u
        transformed = u @ v @ u.conj().T
        off = transformed - torch.diag(torch.diagonal(transformed))
        loss = (off.abs()**2).sum().real/16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    initial_loss = objective(initial.ravel())[0]
    opt = minimize(objective, initial.ravel(), method="L-BFGS-B", jac=True,
                   options={"maxiter": maxiter, "ftol": 1e-16,
                            "gtol": 1e-12, "maxls": 40})
    return opt.x.reshape(14, 4, 3), {"initial_loss": initial_loss,
            "loss": float(opt.fun), "nit": int(opt.nit),
            "nfev": int(opt.nfev), "optimizer_success": bool(opt.success),
            "message": str(opt.message)}


def candidate(edges, angles):
    gates = []
    for slot in range(14):
        for wire in range(4):
            gates.append({"gate": "u3", "qubit": wire,
                          **axis_to_u3(*angles[slot, wire])})
        if slot < 13:
            a, b = edges[slot]
            gates.append({"gate": "cx", "control": a, "target": b})
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=25040)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--warm-sigma", type=float, default=0.05)
    parser.add_argument("--random-sigma", type=float, default=0.4)
    args = parser.parse_args()
    warm = warm_angles()
    rng = np.random.default_rng(args.seed)
    rows = []
    for index, (prefix, tail) in enumerate(zip(PREFIXES, TAILS)):
        ep = endpoint(prefix)
        assert 15 in ep
        assert min(sum(a != b for a, b in zip(prefix,
                       EXACT_PREFIX[:deleted]+EXACT_PREFIX[deleted+1:]))
                   for deleted in range(6)) >= 2
        edges = prefix+tail
        assert len(edges) == 13
        for kind, initial in (("warm", warm+rng.normal(0, args.warm_sigma, warm.shape)),
                              ("random", rng.normal(0, args.random_sigma, warm.shape))):
            angles, record = fit(edges, initial, args.maxiter)
            circuit = candidate(edges, angles)
            path = ROOT / f"search13_joint_gauge_five_endpoint_{index}_{kind}.json"
            path.write_text(json.dumps(circuit, indent=2)+"\n")
            checked = evaluate(circuit)
            record.update(case=index, start=kind, prefix=prefix,
                          endpoint=ep, tail=tail, candidate=path.name,
                          off_diagonal_error=checked["off_diagonal_error"],
                          valid_diagonalizer=checked["valid_diagonalizer"])
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    result = {"warm_source": WARM.name, "seed": args.seed,
              "warm_sigma": args.warm_sigma, "random_sigma": args.random_sigma,
              "maxiter": args.maxiter, "rows": rows,
              "best_loss": min(rows, key=lambda r: r["loss"]),
              "best_off_diagonal": min(rows, key=lambda r: r["off_diagonal_error"])}
    (ROOT / "search13_joint_gauge_five_endpoint_result.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
