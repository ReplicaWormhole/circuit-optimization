"""Read-only KAK angle audit for numerical fourteen-CNOT diagonalizers.

This reports angles after the exact common-Rz center gauge.  Rational-looking
floats are hypotheses, not certificates of an exact circuit.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit.synthesis import OneQubitEulerDecomposer, TwoQubitWeylDecomposition

from topology14_nearhit_kak import factor
from topology14_twolevel_boundary_rank import circuit_matrix
from topology14_kak_center_gauge import rz


EULER = OneQubitEulerDecomposer("ZYZ")


def angle_record(matrix, denominator=16):
    angles = [float(a) for a in EULER.angles(matrix)]
    step = math.pi / denominator
    return {"angles": angles,
            "nearest_numerator": [int(round(a / step)) for a in angles],
            "errors": [float(a - round(a / step) * step) for a in angles]}


def audit(path):
    gates = json.loads(path.read_text())["gates"]
    a02, a13, _ = factor(circuit_matrix(gates[13:65]), (0, 2))
    b01, b23, _ = factor(circuit_matrix(gates[65:]), (0, 1))
    d02, d13, d01, d23 = [TwoQubitWeylDecomposition(u, fidelity=1.0)
                          for u in (a02, a13, b01, b23)]
    n = {0: d02.K2l, 2: d02.K2r, 1: d13.K2l, 3: d13.K2r}
    a_out = {0: d02.K1l, 2: d02.K1r, 1: d13.K1l, 3: d13.K1r}
    b_in = {0: d01.K2l, 1: d01.K2r, 2: d23.K2l, 3: d23.K2r}
    o = {0: d01.K1l, 1: d01.K1r, 2: d23.K1l, 3: d23.K1r}
    center = {q: b_in[q] @ a_out[q] for q in range(4)}
    phi = {q: angle_record(center[q])["angles"][1] for q in range(4)}
    alpha = {(0, 1): phi[0], (2, 3): phi[2]}
    for pair, value in alpha.items():
        for q in pair:
            o[q] = o[q] @ rz(value)
            center[q] = rz(-value) @ center[q]
    frames = {"input": n, "center": center, "output": o}
    records = {name: {str(q): angle_record(m)
                      for q, m in matrices.items()}
               for name, matrices in frames.items()}
    cartan = {str(pair): [float(d.a), float(d.b), float(d.c)]
              for pair, d in zip(((0, 2), (1, 3), (0, 1), (2, 3)),
                                 (d02, d13, d01, d23))}
    return {"source": str(path), "cartan": cartan,
            "center_phase_difference_mod_2pi": {
                "01": float(np.angle(np.exp(1j * (phi[1] - phi[0])))),
                "23": float(np.angle(np.exp(1j * (phi[3] - phi[2]))))},
            "frames": records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = [audit(path) for path in args.input]
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for item in result:
        print(item["source"])
        print("cartan", item["cartan"])
        print("center phase difference", item["center_phase_difference_mod_2pi"])
        for frame, data in item["frames"].items():
            print(frame, [data[str(q)]["nearest_numerator"] for q in range(4)])
            print(frame, [max(map(abs, data[str(q)]["errors"])) for q in range(4)])


if __name__ == "__main__":
    main()
