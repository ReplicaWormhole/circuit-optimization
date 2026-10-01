"""Audit numerical U3 angles in the 16-CNOT Qiskit candidate.

Reports distance to low-denominator rational multiples of pi and recognizes
simple radicals in trigonometric values. This is diagnostic, not a proof that
an angle is irrational or that the circuit is exact.
"""

import argparse
import json
import math
from fractions import Fraction

import sympy as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate")
    parser.add_argument("--max-denominator", type=int, default=64)
    args = parser.parse_args()
    candidate = json.load(open(args.candidate))
    basis = [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5),
             sp.sqrt(2 + sp.sqrt(2)), sp.sqrt(2 - sp.sqrt(2))]
    records = []
    for position, gate in enumerate(candidate["gates"]):
        if gate["gate"] != "u3":
            continue
        for key in ("theta", "phi", "lam"):
            value = gate[key]
            if not isinstance(value, (float, int)):
                continue
            nearest = Fraction(value / math.pi).limit_denominator(args.max_denominator)
            error = abs(value - math.pi * float(nearest))
            if error < 1e-10:
                continue
            record = {"gate_position": position, "qubit": gate["qubit"],
                      "parameter": key, "value": value,
                      "nearest_pi_fraction": str(nearest),
                      "nearest_angle_error": error,
                      "cos": float(math.cos(value)), "sin": float(math.sin(value))}
            for trig in ("cos", "sin"):
                try:
                    guess = sp.nsimplify(record[trig], basis, tolerance=1e-13,
                                         full=False)
                    if guess.is_number and abs(float(guess) - record[trig]) < 1e-12:
                        record[trig + "_radical_guess"] = str(guess)
                except (ValueError, TypeError):
                    pass
            records.append(record)
    print(json.dumps({"candidate": args.candidate,
                      "non_pi_fraction_angles_up_to_denominator": args.max_denominator,
                      "count": len(records), "angles": records}, indent=2))


if __name__ == "__main__":
    main()
