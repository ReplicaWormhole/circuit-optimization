"""Finite-difference Gauss-Newton refinement of a near-exact 14-CNOT list.

The output eigenvalue labels are unambiguous at the supplied near-hit.  We
hold its six-CNOT exact prefix and eight-CNOT topology fixed, then solve the
matrix equation U V4 = D U in the 108 local ZYZ angles.  This is numerical
evidence, not exact angle certification.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift


ROOT = Path(__file__).resolve().parent
ROOTS = np.array([1, 1j, -1, -1j])
V = right_shift(4)


def parse(source):
    gates = source["gates"]
    prefix = gates[:13]
    suffix = gates[13:]
    if sum(gate["gate"] == "cx" for gate in prefix) != 6:
        raise ValueError("expected exact six-CNOT prefix")
    if len(suffix) != 9 * 12 + 8:
        raise ValueError(f"unexpected suffix gate count {len(suffix)}")
    angles = np.zeros((9, 4, 3), dtype=float)
    edges = []
    index = 0
    for slot in range(9):
        for qubit in range(4):
            for axis, name in enumerate(("rz", "ry", "rz")):
                gate = suffix[index]
                if gate["gate"] != name or gate["qubit"] != qubit:
                    raise ValueError(f"unexpected layer gate {index}: {gate}")
                angles[slot, qubit, axis] = float(gate["theta"])
                index += 1
        if slot < 8:
            gate = suffix[index]
            if gate["gate"] != "cx":
                raise ValueError(f"expected CNOT at suffix index {index}")
            edges.append((gate["control"], gate["target"]))
            index += 1
    prefix_u = np.eye(16, dtype=complex)
    for gate in prefix:
        matrix = (cnot(gate["control"], gate["target"], 4)
                  if gate["gate"] == "cx" else
                  embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        prefix_u = matrix @ prefix_u
    return prefix, prefix_u, edges, angles


def rz(angle):
    return np.diag(np.exp(np.array([-0.5j, 0.5j]) * angle))


def ry(angle):
    c, s = math.cos(angle / 2), math.sin(angle / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def unitary(prefix_u, cx_matrices, angles):
    out = prefix_u
    for slot in range(9):
        wires = [rz(angles[slot, q, 2]) @ ry(angles[slot, q, 1])
                 @ rz(angles[slot, q, 0]) for q in range(4)]
        layer = wires[0]
        for wire in wires[1:]:
            layer = np.kron(layer, wire)
        out = layer @ out
        if slot < 8:
            out = cx_matrices[slot] @ out
    return out


def residual(u, labels):
    difference = u @ V - labels[:, None] * u
    return np.concatenate((difference.real.ravel(), difference.imag.ravel()))


def emit(prefix, edges, angles):
    gates = list(prefix)
    for slot in range(9):
        for qubit in range(4):
            for axis, name in enumerate(("rz", "ry", "rz")):
                gates.append({"gate": name, "qubit": qubit,
                              "theta": float(angles[slot, qubit, axis])})
        if slot < 8:
            gates.append({"gate": "cx", "control": edges[slot][0],
                          "target": edges[slot][1]})
    return {"n": 4, "gates": gates}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    p.add_argument("--maxiter", type=int, default=6)
    p.add_argument("--difference-step", type=float, default=1e-5)
    p.add_argument("--rcond", type=float, default=1e-10)
    args = p.parse_args()
    source = json.loads(args.input.read_text())
    prefix, prefix_u, edges, angles = parse(source)
    cx_matrices = [cnot(*edge, 4) for edge in edges]
    first_u = unitary(prefix_u, cx_matrices, angles)
    diagonal = np.diag(first_u @ V @ first_u.conj().T)
    root_indices = np.argmin(abs(diagonal[:, None] - ROOTS[None, :]), axis=1)
    labels = ROOTS[root_indices]
    if np.max(abs(diagonal - labels)) > 1e-5:
        raise ValueError("near-hit output labels are not unambiguous")
    history = []
    for iteration in range(args.maxiter):
        r = residual(unitary(prefix_u, cx_matrices, angles), labels)
        old_norm = float(np.linalg.norm(r))
        old_max = float(np.max(abs(r)))
        if old_max < 1e-13:
            history.append({"iteration": iteration, "residual_norm": old_norm,
                            "max_component": old_max, "stopped": "threshold"})
            break
        flat = angles.ravel()
        jacobian = np.empty((len(r), len(flat)), dtype=float)
        for coordinate in range(len(flat)):
            plus, minus = flat.copy(), flat.copy()
            plus[coordinate] += args.difference_step
            minus[coordinate] -= args.difference_step
            rp = residual(unitary(prefix_u, cx_matrices,
                                  plus.reshape(9, 4, 3)), labels)
            rm = residual(unitary(prefix_u, cx_matrices,
                                  minus.reshape(9, 4, 3)), labels)
            jacobian[:, coordinate] = (rp - rm) / (2 * args.difference_step)
        delta, _, rank, singular = np.linalg.lstsq(
            jacobian, -r, rcond=args.rcond)
        if np.linalg.norm(delta) > 0.1:
            delta *= 0.1 / np.linalg.norm(delta)
        accepted = False
        for factor in (1.0, 0.5, 0.25, 0.125):
            trial = (flat + factor * delta).reshape(9, 4, 3)
            trial_norm = float(np.linalg.norm(residual(
                unitary(prefix_u, cx_matrices, trial), labels)))
            if trial_norm < old_norm:
                angles = trial
                accepted = True
                break
        history.append({"iteration": iteration, "residual_norm": old_norm,
                        "max_component": old_max, "jacobian_rank": int(rank),
                        "smallest_retained_singular": float(singular[rank - 1]),
                        "step_norm": float(np.linalg.norm(delta)),
                        "accepted_factor": factor if accepted else None,
                        "next_norm": trial_norm if accepted else None})
        if not accepted:
            break
    candidate = emit(prefix, edges, angles)
    args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    final_r = residual(unitary(prefix_u, cx_matrices, angles), labels)
    result = {"source": str(args.input), "maxiter": args.maxiter,
              "difference_step": args.difference_step, "rcond": args.rcond,
              "output_labels": root_indices.tolist(),
              "initial_offdiag": evaluate(source)["off_diagonal_error"],
              "final_residual_norm": float(np.linalg.norm(final_r)),
              "final_max_component": float(np.max(abs(final_r))),
              "final_offdiag": check["off_diagonal_error"],
              "valid_diagonalizer": check["valid_diagonalizer"],
              "cnot_count": check["cnot_count"], "history": history}
    args.result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
