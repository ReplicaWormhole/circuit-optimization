"""Continue the checker-valid p≈3/4 mixer toward p=1/2 in small steps."""

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
TARGETS = (0.70, 0.65, 0.60, 0.55, 0.50)
torch.set_num_threads(1)


def main():
    source = json.loads((ROOT / "algebraic14_mixer_p75.json").read_text())
    old = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix_gates = source["gates"][:13]
    initial, edges = layers_from_suffix(source["gates"][13:])
    prefix = torch.tensor(Operator(to_qiskit(
        {"n": 4, "gates": prefix_gates})).data, dtype=DTYPE)
    old_u = torch.tensor(Operator(to_qiskit(old)).data, dtype=DTYPE)
    v = torch.tensor(right_shift(4), dtype=DTYPE)
    labels = json.loads(
        (ROOT / "algebraic14_pure_orbit_refine_result.json").read_text())["labels"]
    roots = np.array([1, 1j, -1, -1j])
    d = torch.diag(torch.tensor(roots[labels], dtype=DTYPE))
    stages = []
    for target in TARGETS:
        def objective(flat):
            parameters = torch.tensor(flat.reshape(initial.shape),
                                      dtype=torch.float64, requires_grad=True)
            u = circuit(parameters, edges, prefix)
            residual = u @ v - d @ u
            diag_loss = (torch.abs(residual) ** 2).sum().real / 16
            m = u @ old_u.conj().T
            p1 = torch.abs(m[5, 6]) ** 2
            p2 = torch.abs(m[6, 5]) ** 2
            gauge_loss = (p1 - target) ** 2 + (p2 - target) ** 2
            total = 1e8 * diag_loss + gauge_loss
            (gradient,) = torch.autograd.grad(total, parameters)
            return float(total.detach()), gradient.detach().numpy().ravel().copy()

        optimized = minimize(objective, initial.ravel(),
                             method="L-BFGS-B", jac=True,
                             options={"maxiter": 500, "ftol": 1e-17,
                                      "gtol": 1e-11, "maxls": 40})
        initial = optimized.x.reshape(initial.shape)
        candidate = emit(prefix_gates, initial, edges)
        path = ROOT / f"algebraic14_mixer_cont_p{int(target*100):02d}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        u = Operator(to_qiskit(candidate)).data
        old_np = Operator(to_qiskit(old)).data
        m = u @ old_np.conj().T
        stages.append({"target_p": target,
                       "achieved_p": float(abs(m[5, 6]) ** 2),
                       "offdiag_error": check["off_diagonal_error"],
                       "valid": check["valid_diagonalizer"],
                       "iterations": int(optimized.nit),
                       "objective": float(optimized.fun),
                       "candidate": str(path)})
    result = {"stages": stages}
    (ROOT / "algebraic14_mixer_continuation_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
