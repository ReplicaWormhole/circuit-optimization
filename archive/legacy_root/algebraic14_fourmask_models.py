"""Enumerate four-mask phase-polynomial supports modulo shift-invariant gauge.

Every coefficient is an integer multiple of pi/8 modulo 2pi. For each
four-element nonlocal support, solve the 16 orbit-edge phase equations with
arbitrary local Z terms. Write one model per support for subsequent parity
network searches. This is finite exact modular arithmetic via Z3.
"""

import json
from itertools import combinations
from pathlib import Path

from z3 import Int, Solver, sat, unsat

from algebraic14_phase_support import NONLOCAL, SINGLES, TARGET, signature


ROOT = Path(__file__).parent


def model_for_support(support):
    masks = SINGLES + support
    coeff = {mask: Int(f"a_{mask}") for mask in masks}
    solver = Solver()
    solver.set("timeout", 10000)
    for variable in coeff.values():
        solver.add(variable >= 0, variable < 16)
    for mask in support:
        solver.add(coeff[mask] != 0)
    for x in range(16):
        lhs = sum(coeff[mask] * signature(mask, x) for mask in masks)
        rhs = sum(value * signature(mask, x)
                  for mask, value in TARGET.items())
        solver.add((lhs - rhs) % 16 == 0)
    status = solver.check()
    if status == sat:
        assignment = solver.model()
        return {str(mask): assignment[variable].as_long()
                for mask, variable in coeff.items()
                if assignment[variable].as_long() != 0}
    if status != unsat:
        raise RuntimeError(f"solver returned {status} for support {support}")
    return None


def main():
    models = []
    for support in combinations(NONLOCAL, 4):
        coefficients = model_for_support(support)
        if coefficients is not None:
            models.append({"support": support, "coefficients_pi_over_8": coefficients})
    path = ROOT / "algebraic14_fourmask_models.json"
    path.write_text(json.dumps(models, indent=2) + "\n")
    print(json.dumps({"tested_supports": 210, "satisfiable_supports": len(models),
                      "models_file": path.name,
                      "supports": [record["support"] for record in models]}, indent=2))


if __name__ == "__main__":
    main()
