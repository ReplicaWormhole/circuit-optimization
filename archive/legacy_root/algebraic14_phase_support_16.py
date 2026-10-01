"""Check whether a pi/16 angle grid changes minimum parity support.

The baseline B4*CZ phase coefficients are doubled from pi/8 units to pi/16
units. Z3 tests whether <=3 nonlocal parity masks can reproduce the orbit
phase differences modulo 2pi, with arbitrary local terms. A 4-mask SAT
control guards against a convention mistake. This is one bounded grid test.
"""

import json

from z3 import AtMost, Int, Solver, sat, unsat

from algebraic14_phase_support import NONLOCAL, SINGLES, TARGET, signature


def solve(max_nonlocal):
    masks = SINGLES + NONLOCAL
    coeff = {mask: Int(f"a_{mask}") for mask in masks}
    solver = Solver()
    solver.set("timeout", 30000)
    for variable in coeff.values():
        solver.add(variable >= 0, variable < 32)
    solver.add(AtMost(*[coeff[mask] != 0 for mask in NONLOCAL], max_nonlocal))
    for x in range(16):
        lhs = sum(coeff[mask] * signature(mask, x) for mask in masks)
        rhs = sum(2 * value * signature(mask, x)
                  for mask, value in TARGET.items())
        solver.add((lhs - rhs) % 32 == 0)
    status = solver.check()
    model = None
    if status == sat:
        answer = solver.model()
        model = {str(mask): answer[variable].as_long()
                 for mask, variable in coeff.items()
                 if answer[variable].as_long()}
    return {"max_nonlocal": max_nonlocal,
            "status": str(status), "model": model,
            "reason_unknown": solver.reason_unknown()
            if status not in (sat, unsat) else None}


def main():
    print(json.dumps({"angle_grid": "pi/16", "results": [solve(3), solve(4)]},
                     indent=2))


if __name__ == "__main__":
    main()
