"""Refine the 14-CX p=1/2 output mixer without a gauge penalty.

The prior constrained run reached p≈1/2 with V4 error 9e-9. Starting
there, minimize only fixed-label UV4−DU residual, then measure whether
p remains close to 1/2. This tests if the simple mixer lies on the exact
diagonalizer manifold; it does not certify exact gates.
"""

import json
from pathlib import Path

import numpy as np
import torch
from qiskit.quantum_info import Operator
from scipy.optimize import minimize

from check_circuit import evaluate, right_shift
from qiskit_crosscheck import to_qiskit
from topology14_eightcx_opt import emit, layers_from_suffix
from topology14_fourcx_opt import DTYPE, circuit


ROOT = Path(__file__).parent
torch.set_num_threads(1)


def main():
    source = json.loads((ROOT / "algebraic14_mixer_p50.json").read_text())
    old = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix_gates = source["gates"][:13]
    initial, edges = layers_from_suffix(source["gates"][13:])
    prefix = torch.tensor(Operator(to_qiskit(
        {"n": 4, "gates": prefix_gates})).data, dtype=DTYPE)
    v = torch.tensor(right_shift(4), dtype=DTYPE)
    labels = json.loads(
        (ROOT / "algebraic14_pure_orbit_refine_result.json").read_text())["labels"]
    roots = np.array([1, 1j, -1, -1j])
    d = torch.diag(torch.tensor(roots[labels], dtype=DTYPE))

    def objective(flat):
        parameters = torch.tensor(flat.reshape(initial.shape),
                                  dtype=torch.float64, requires_grad=True)
        u = circuit(parameters, edges, prefix)
        residual = u @ v - d @ u
        loss = (torch.abs(residual) ** 2).sum().real / 16
        (gradient,) = torch.autograd.grad(loss, parameters)
        return float(loss.detach() * 1e12), (
            gradient.detach().numpy().ravel().copy() * 1e12)

    optimized = minimize(objective, initial.ravel(),
                         method="L-BFGS-B", jac=True,
                         options={"maxiter": 1000, "ftol": 1e-17,
                                  "gtol": 1e-11, "maxls": 40})
    candidate = emit(prefix_gates,
                     optimized.x.reshape(initial.shape), edges)
    path = ROOT / "algebraic14_mixer_p50_refined.json"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    u = Operator(to_qiskit(candidate)).data
    u15 = Operator(to_qiskit(old)).data
    m = u @ u15.conj().T
    result = {"initial_scaled_loss": objective(initial.ravel())[0],
              "final_scaled_loss": float(optimized.fun),
              "iterations": int(optimized.nit),
              "cnot_count": check["cnot_count"],
              "offdiag_error": check["off_diagonal_error"],
              "valid": check["valid_diagonalizer"],
              "mix_probabilities": [float(abs(m[5, 6]) ** 2),
                                    float(abs(m[6, 5]) ** 2)],
              "candidate": str(path)}
    (ROOT / "algebraic14_mixer_p50_refine_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
