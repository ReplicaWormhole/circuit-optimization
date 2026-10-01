"""Test whether 14-CX degenerate two-level mixer amplitudes are freely tunable.

Starting from run148's exact numerical diagonalizer, optimize the same
six-prefix/eight-tail topology with fixed V4 eigenvalue labels and a soft
constraint |(U U15†)[5,6]|² = p for p=1/4,1/2,3/4. The paired [6,5]
entry is constrained too. If both V4 residual and gauge penalty vanish,
the mixer can be chosen at a simple angle. This is numerical evidence
only; exact local gate specification remains separate.
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
TARGETS = (0.25, 0.5, 0.75)
torch.set_num_threads(1)


def main():
    source = json.loads(
        (ROOT / "algebraic14_pure_orbit_refined_candidate.json").read_text())
    exact15 = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    prefix_gates = source["gates"][:13]
    initial, edges = layers_from_suffix(source["gates"][13:])
    prefix = torch.tensor(Operator(to_qiskit(
        {"n": 4, "gates": prefix_gates})).data, dtype=DTYPE)
    old = torch.tensor(Operator(to_qiskit(exact15)).data, dtype=DTYPE)
    v = torch.tensor(right_shift(4), dtype=DTYPE)
    labels = json.loads(
        (ROOT / "algebraic14_pure_orbit_refine_result.json").read_text())["labels"]
    roots = np.array([1, 1j, -1, -1j])
    d = torch.diag(torch.tensor(roots[labels], dtype=DTYPE))
    results = []
    for target in TARGETS:
        calls = [0]

        def objective(flat):
            parameters = torch.tensor(flat.reshape(initial.shape),
                                      dtype=torch.float64,
                                      requires_grad=True)
            u = circuit(parameters, edges, prefix)
            residual = u @ v - d @ u
            diag_loss = (torch.abs(residual) ** 2).sum().real / 16
            mixing = u @ old.conj().T
            amplitudes = torch.stack(
                (torch.abs(mixing[5, 6]) ** 2,
                 torch.abs(mixing[6, 5]) ** 2))
            gauge_loss = torch.sum((amplitudes - target) ** 2)
            total = 1e8 * diag_loss + gauge_loss
            (gradient,) = torch.autograd.grad(total, parameters)
            calls[0] += 1
            return float(total.detach()), (
                gradient.detach().numpy().ravel().copy())

        optimized = minimize(objective, initial.ravel(),
                             method="L-BFGS-B", jac=True,
                             options={"maxiter": 600, "ftol": 1e-16,
                                      "gtol": 1e-11, "maxls": 40})
        candidate = emit(prefix_gates,
                         optimized.x.reshape(initial.shape), edges)
        path = ROOT / f"algebraic14_mixer_p{int(target * 100):02d}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        u = Operator(to_qiskit(candidate)).data
        u15 = Operator(to_qiskit(exact15)).data
        m = u @ u15.conj().T
        results.append({
            "target_mix_probability": target,
            "achieved_mix_probability": [
                float(abs(m[5, 6]) ** 2),
                float(abs(m[6, 5]) ** 2)],
            "cnot_count": check["cnot_count"],
            "offdiag_error": check["off_diagonal_error"],
            "valid": check["valid_diagonalizer"],
            "iterations": int(optimized.nit),
            "evaluations": calls[0],
            "objective": float(optimized.fun),
            "candidate": str(path),
        })
    result = {"trials": results}
    (ROOT / "algebraic14_mixer_gauge_scan_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
