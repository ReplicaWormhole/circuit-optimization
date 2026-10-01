"""Ask whether B4*CZ admits a <=3-mask nonlocal parity phase polynomial.

Angles are integral multiples of pi/8, modulo 2pi; unrestricted single-wire
Z rotations are allowed on the same angle grid. Equality is tested only up
to multiplication by a shift-invariant diagonal unitary, by comparing each
phase difference along the V4 orbit. This is a bounded algebraic ansatz, not
a circuit lower bound or a statement about arbitrary real angles.
"""

import json
from z3 import AtMost, Int, Solver, sat, unsat


SINGLES = (8, 4, 2, 1)
NONLOCAL = tuple(mask for mask in range(1, 16)
                 if mask.bit_count() >= 2 and mask != 15)
TARGET = {4: 5, 2: 6, 1: 3, 11: 1, 13: 2, 14: 3, 6: -4}


def shift(x):
    return ((x & 1) << 3) | (x >> 1)


def signature(mask, x):
    parity = lambda y: (-1) ** ((mask & y).bit_count() % 2)
    return (parity(shift(x)) - parity(x)) // 2


def solve(max_nonlocal):
    coeff = {mask: Int(f"a_{mask}") for mask in SINGLES + NONLOCAL}
    solver = Solver()
    solver.set("timeout", 30000)
    for variable in coeff.values():
        solver.add(variable >= 0, variable < 16)
    solver.add(AtMost(*[coeff[mask] != 0 for mask in NONLOCAL], max_nonlocal))
    for x in range(16):
        left = sum(coeff[mask] * signature(mask, x) for mask in coeff)
        right = sum(value * signature(mask, x)
                    for mask, value in TARGET.items())
        solver.add((left - right) % 16 == 0)
    status = solver.check()
    model = None
    if status == sat:
        answer = solver.model()
        model = {str(mask): answer[variable].as_long()
                 for mask, variable in coeff.items()
                 if answer[variable].as_long()}
    return {"max_nonlocal": max_nonlocal, "status": str(status),
            "model": model, "reason_unknown": solver.reason_unknown()
            if status not in (sat, unsat) else None}


def main():
    results = [solve(3), solve(4)]
    print(json.dumps({"angle_grid": "pi/8", "nonlocal_masks": NONLOCAL,
                      "results": results}, indent=2))


if __name__ == "__main__":
    main()
