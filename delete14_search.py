"""Bounded CNOT deletion plus edge mutation from the exact 15-CNOT circuit.

The label-free objective is the off-diagonal Frobenius norm of U V U†.
Corrections are arbitrary SU(2) axis-angle rotations on every wire between
successive CNOTs. This is a local numerical search, not an exact certificate.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift


ROOT = Path(__file__).resolve().parent
torch.set_num_threads(1)
CDTYPE = torch.complex128
RDTYPE = torch.float64


def fixed_matrix(gate):
    if gate["gate"] == "cx":
        return cnot(gate["control"], gate["target"], 4)
    return embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4)


def split_baseline(gates, deleted_cx):
    chunks = [[]]
    cx_gates = []
    cx_idx = -1
    for gate in gates:
        if gate["gate"] == "cx":
            cx_idx += 1
            if cx_idx == deleted_cx:
                continue
            cx_gates.append(gate)
            chunks.append([])
        else:
            chunks[-1].append(gate)
    if cx_idx != 14:
        raise ValueError(f"expected 15 CNOTs, found {cx_idx + 1}")
    return chunks, cx_gates


def kron_gate(mat, qubit):
    result = torch.ones((1, 1), dtype=CDTYPE)
    identity = torch.eye(2, dtype=CDTYPE)
    for j in range(4):
        result = torch.kron(result, mat if j == qubit else identity)
    return result


def axis_rotation(x, y, z):
    # exp[-i(x X + y Y + z Z)/2], regular at the identity.
    # The tiny regularizer avoids the undefined AD derivative of sqrt at zero.
    length = torch.sqrt(x*x + y*y + z*z + 1e-24)
    c = torch.cos(length / 2)
    s_over_length = 0.5 * torch.sinc(length / (2 * math.pi))
    return torch.stack((torch.stack((c - 1j*z*s_over_length,
                                     (-1j*x-y)*s_over_length)),
                        torch.stack(((-1j*x+y)*s_over_length,
                                     c + 1j*z*s_over_length))))


def axis_to_u3(x, y, z):
    length = math.sqrt(x*x + y*y + z*z)
    c = math.cos(length / 2)
    s_over_length = 0.5 * np.sinc(length / (2 * math.pi))
    a = c - 1j*z*s_over_length
    lower = (-1j*x+y)*s_over_length
    upper = (-1j*x-y)*s_over_length
    theta = 2 * math.atan2(abs(lower), abs(a))
    alpha = np.angle(a)
    if abs(lower) < 1e-14:
        # At theta=0, the off-diagonal phases are undefined.  Choose a
        # diagonal U3 matching the SU(2) matrix up to one global phase.
        return {"theta": 0.0, "phi": 0.0, "lam": float(-2 * alpha)}
    phi = np.angle(lower) - alpha
    lam = np.angle(-upper) - alpha
    return {"theta": theta, "phi": float(phi), "lam": float(lam)}


def u3_to_axis(gate):
    theta, phi, lam = (float(gate[k]) for k in ("theta", "phi", "lam"))
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    a = c * np.exp(-0.5j * (phi + lam))
    b = -s * np.exp(0.5j * (lam - phi))
    components = np.array([-b.imag, -b.real, -a.imag], dtype=float)
    sine = np.linalg.norm(components)
    if sine < 1e-14:
        return np.zeros(3)
    length = 2 * math.atan2(sine, a.real)
    return components * (length / sine)


def setup(deleted_cx, base_candidate, replace_cx=None, replacement=None):
    base = json.loads(base_candidate.read_text())
    chunks, cxs = split_baseline(base["gates"], deleted_cx)
    if replace_cx is not None:
        if not 0 <= replace_cx < len(cxs):
            raise ValueError("replacement CNOT index out of range")
        if replacement is None or replacement[0] == replacement[1]:
            raise ValueError("replacement must specify distinct wires")
        cxs[replace_cx] = {"gate": "cx", "control": replacement[0],
                           "target": replacement[1]}
    matrices = []
    for chunk in chunks:
        mat = np.eye(16, dtype=complex)
        for gate in chunk:
            mat = fixed_matrix(gate) @ mat
        matrices.append(torch.as_tensor(mat, dtype=CDTYPE))
    cx_matrices = [torch.as_tensor(fixed_matrix(g), dtype=CDTYPE) for g in cxs]
    vb = right_shift(4)
    baseline_u = np.eye(16, dtype=complex)
    for gate in base["gates"]:
        baseline_u = fixed_matrix(gate) @ baseline_u
    diag = np.diag(baseline_u @ vb @ baseline_u.conj().T)
    roots = np.array([1, 1j, -1, -1j])
    exact_diag = roots[np.argmin(np.abs(diag[:, None] - roots[None, :]), axis=1)]
    return chunks, cxs, matrices, cx_matrices, torch.as_tensor(vb, dtype=CDTYPE), torch.diag(torch.as_tensor(exact_diag, dtype=CDTYPE))


def optimize_one(deleted_cx, seed, maxiter, output, initial_candidate=None,
                 optimizer="L-BFGS-B", objective_mode="fixed-label", base_candidate=None,
                 initial_sigma=0.01, replace_cx=None, replacement=None):
    if base_candidate is None:
        base_candidate = ROOT / "topology16_15_exact_candidate.json"
    chunks, cxs, mats, cx_mats, v, d = setup(deleted_cx, base_candidate,
                                            replace_cx, replacement)
    nslots = len(chunks)
    rng = np.random.default_rng(seed)
    x0 = np.zeros((nslots, 4, 3), dtype=np.float64)
    if initial_candidate:
        old = json.loads(initial_candidate.read_text())["gates"]
        position = 0
        for slot, chunk in enumerate(chunks):
            if old[position:position + len(chunk)] != chunk:
                raise ValueError("warm-start candidate has different fixed gates")
            position += len(chunk)
            for q in range(4):
                gate = old[position]
                if gate["gate"] != "u3" or gate["qubit"] != q:
                    raise ValueError("warm-start correction layout differs")
                x0[slot, q] = u3_to_axis(gate)
                position += 1
            if slot < len(cxs):
                if old[position] != cxs[slot]:
                    raise ValueError("warm-start CNOT topology differs")
                position += 1
        if position != len(old):
            raise ValueError("warm-start candidate has extra gates")
    elif seed:
        x0 += rng.normal(0, initial_sigma, x0.shape)
    history = []

    def objective(flat):
        angles = torch.tensor(flat.reshape(nslots, 4, 3), dtype=RDTYPE, requires_grad=True)
        u = torch.eye(16, dtype=CDTYPE)
        for slot in range(nslots):
            u = mats[slot] @ u
            for q in range(4):
                a, b, c = angles[slot, q]
                u = kron_gate(axis_rotation(a, b, c), q) @ u
            if slot < len(cx_mats):
                u = cx_mats[slot] @ u
        if objective_mode == "fixed-label":
            diff = u @ v - d @ u
            loss = (diff.abs() ** 2).sum().real / 16
        else:
            transformed = u @ v @ u.conj().T
            offdiag = transformed - torch.diag(torch.diagonal(transformed))
            loss = (offdiag.abs() ** 2).sum().real / 16
        grad = torch.autograd.grad(loss, angles)[0]
        return float(loss.detach()), grad.detach().numpy().ravel().copy()

    def callback(x):
        if len(history) % 25 == 0:
            history.append(objective(x)[0])
        else:
            history.append(None)

    initial = objective(x0.ravel())[0]
    options = ({"maxiter": maxiter, "ftol": 1e-16, "gtol": 1e-13, "maxls": 30}
               if optimizer == "L-BFGS-B" else
               {"maxiter": maxiter, "gtol": 1e-13})
    result = minimize(objective, x0.ravel(), method=optimizer, jac=True,
                      options=options, callback=callback)
    angles = result.x.reshape(nslots, 4, 3)
    gates = []
    for slot, chunk in enumerate(chunks):
        gates.extend(chunk)
        for q in range(4):
            a, b, c = angles[slot, q]
            gates.append({"gate": "u3", "qubit": q, **axis_to_u3(a, b, c)})
        if slot < len(cxs):
            gates.append(cxs[slot])
    candidate = {"n": 4, "gates": gates}
    output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    print(json.dumps({"deleted_cx": deleted_cx, "replace_cx": replace_cx,
                      "replacement": replacement,
                      "seed": seed,
                      "initial_sigma": initial_sigma,
                      "maxiter": maxiter, "nit": result.nit, "nfev": result.nfev,
                      "initial_loss": initial, "final_loss": result.fun,
                      "optimizer": optimizer,
                      "objective_mode": objective_mode,
                      "optimizer_success": bool(result.success),
                      "optimizer_message": str(result.message),
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid_diagonalizer": check["valid_diagonalizer"],
                      "candidate": str(output)}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--deleted-cx", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--initial-sigma", type=float, default=0.01)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--initial-candidate", type=Path)
    parser.add_argument("--optimizer", choices=["L-BFGS-B", "BFGS"], default="L-BFGS-B")
    parser.add_argument("--objective", choices=["fixed-label", "off-diagonal"],
                        default="fixed-label")
    parser.add_argument("--base-candidate", type=Path,
                        default=ROOT / "topology16_15_exact_candidate.json")
    parser.add_argument("--replace-cx", type=int)
    parser.add_argument("--replacement-control", type=int)
    parser.add_argument("--replacement-target", type=int)
    args = parser.parse_args()
    optimize_one(args.deleted_cx, args.seed, args.maxiter, args.output,
                 args.initial_candidate, args.optimizer, args.objective,
                 args.base_candidate, args.initial_sigma, args.replace_cx,
                 (args.replacement_control, args.replacement_target)
                 if args.replace_cx is not None else None)


if __name__ == "__main__":
    main()
