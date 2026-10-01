"""Inspect local KAK boundary gauges of numerical 14-CNOT paired blocks.

For each of two disjoint-pair layers, take Qiskit's 2-qubit Weyl factors.
Combine the four local factors at the boundary between the two Cartan
layers, and compare all one-qubit ZYZ angles with a pi/8 lattice.
"""

import argparse
import json
import math
from pathlib import Path

from qiskit.synthesis import OneQubitEulerDecomposer, TwoQubitWeylDecomposition

from topology14_nearhit_kak import factor
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
DECOMPOSER = OneQubitEulerDecomposer("ZYZ")


def angle_record(matrix):
    theta, phi, lam = [float(value) for value in DECOMPOSER.angles(matrix)]
    angles = [theta, phi, lam]
    ratios = [value / math.pi for value in angles]
    rounded = [round(8 * value / math.pi) for value in angles]
    errors = [float(value - k * math.pi / 8)
              for value, k in zip(angles, rounded)]
    return {"zyz_radians": angles, "zyz_over_pi": ratios,
            "nearest_pi_over_8": rounded,
            "nearest_pi_over_8_errors": errors,
            "max_abs_grid_error": max(map(abs, errors))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    gates = json.loads(args.candidate.read_text())["gates"]
    a02, a13, _ = factor(circuit_matrix(gates[13:65]), (0, 2))
    b01, b23, _ = factor(circuit_matrix(gates[65:]), (0, 1))
    da02, da13, db01, db23 = map(TwoQubitWeylDecomposition,
                                  (a02, a13, b01, b23))
    a_input = {0: da02.K2l, 2: da02.K2r,
               1: da13.K2l, 3: da13.K2r}
    a_output = {0: da02.K1l, 2: da02.K1r,
                1: da13.K1l, 3: da13.K1r}
    b_input = {0: db01.K2l, 1: db01.K2r,
               2: db23.K2l, 3: db23.K2r}
    b_output = {0: db01.K1l, 1: db01.K1r,
                2: db23.K1l, 3: db23.K1r}
    result = {"candidate": str(args.candidate),
              "input_local": {str(q): angle_record(a_input[q]) for q in range(4)},
              "center_local": {str(q): angle_record(b_input[q] @ a_output[q])
                               for q in range(4)},
              "output_local": {str(q): angle_record(b_output[q]) for q in range(4)}}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"input_max_errors": {q: row["max_abs_grid_error"]
                                            for q, row in result["input_local"].items()},
                      "center": result["center_local"],
                      "output_max_errors": {q: row["max_abs_grid_error"]
                                             for q, row in result["output_local"].items()}},
                     indent=2))


if __name__ == "__main__":
    main()
