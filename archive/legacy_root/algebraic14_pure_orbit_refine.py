"""Refine run146's 14-CX near-hit with fixed eigenvalue labels.

The eight-CX tail topology and exact six-CX prefix stay fixed. Reoptimize
all local ZYZ angles by L-BFGS-B against ||U V4 - D U||^2, where D is inferred
from the near-hit. Report full checker/Qiskit numerical residuals. Exact
rational-angle certification is a separate step.
"""

import json
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize
from qiskit.quantum_info import Operator

from check_circuit import evaluate, right_shift
from qiskit_crosscheck import to_qiskit
from topology14_eightcx_opt import emit, layers_from_suffix
from topology14_fourcx_opt import DTYPE, circuit


ROOT = Path(__file__).parent
torch.set_num_threads(1)


def main():
    near = json.loads((ROOT / "algebraic14_pure_orbit_tail8_best.json").read_text())
    prefix_gates = near["gates"][:13]
    suffix = near["gates"][13:]
    initial, edges = layers_from_suffix(suffix)
    assert len(edges) == 8
    prefix_np = Operator(to_qiskit(
        {"n": 4, "gates": prefix_gates})).data
    u_near = Operator(to_qiskit(near)).data
    v_np = right_shift(4)
    diagonal = np.diag(u_near @ v_np @ u_near.conj().T)
    roots = np.array([1, 1j, -1, -1j])
    labels = np.argmin(abs(diagonal[:, None] - roots[None, :]), axis=1)
    d_np = np.diag(roots[labels])
    prefix = torch.tensor(prefix_np, dtype=DTYPE)
    v = torch.tensor(v_np, dtype=DTYPE)
    d = torch.tensor(d_np, dtype=DTYPE)
    evaluations = [0]

    def objective(flat):
        parameters = torch.tensor(flat.reshape(initial.shape),
                                  dtype=torch.float64, requires_grad=True)
        u = circuit(parameters, edges, prefix)
        difference = u @ v - d @ u
        loss = (torch.abs(difference) ** 2).sum().real / 16
        (gradient,) = torch.autograd.grad(loss, parameters)
        evaluations[0] += 1
        return float(loss.detach() * 1e12), (
            gradient.detach().numpy().ravel().copy() * 1e12)

    initial_scaled, _ = objective(initial.ravel())
    optimized = minimize(objective, initial.ravel(),
                         method="L-BFGS-B", jac=True,
                         options={"maxiter": 1000, "ftol": 1e-16,
                                  "gtol": 1e-10, "maxls": 40})
    final_scaled, _ = objective(optimized.x)
    candidate = emit(prefix_gates,
                     optimized.x.reshape(initial.shape), edges)
    path = ROOT / "algebraic14_pure_orbit_refined_candidate.json"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    u_final = Operator(to_qiskit(candidate)).data
    residual = u_final @ v_np - d_np @ u_final
    result = {
        "initial_scaled_loss": initial_scaled,
        "final_scaled_loss": final_scaled,
        "scaling": 1e12,
        "iterations": int(optimized.nit),
        "evaluations": evaluations[0],
        "optimizer_success": bool(optimized.success),
        "optimizer_message": str(optimized.message),
        "labels": [int(label) for label in labels],
        "cnot_count": check["cnot_count"],
        "checker_offdiag": check["off_diagonal_error"],
        "checker_valid": check["valid_diagonalizer"],
        "qiskit_fixed_label_max_residual": float(
            np.max(np.abs(residual))),
        "candidate": str(path),
    }
    (ROOT / "algebraic14_pure_orbit_refine_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
