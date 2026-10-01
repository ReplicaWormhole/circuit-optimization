"""Gauge-invariant two-qubit Weyl coordinates of run146's 14-CX near hit.

The eight-CX tail has two disjoint-pair layers: first pairs (0,2),(1,3),
then (0,1),(2,3). Factor each exact numerical layer across that bipartition
and compute the four two-qubit Cartan triples. These invariants discard the
many unconstrained local ZYZ angles and can reveal rational-pi structure.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit.synthesis import TwoQubitWeylDecomposition

from topology14_boundary_schmidt import bits
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent


def realign(u, left):
    right = tuple(q for q in range(4) if q not in left)
    out = np.zeros((16, 16), dtype=complex)
    for output in range(16):
        for inp in range(16):
            row = bits(output, left) * 4 + bits(inp, left)
            col = bits(output, right) * 4 + bits(inp, right)
            out[row, col] = u[output, inp]
    return out


def factor(u, left):
    realigned = realign(u, left)
    v, singular, vh = np.linalg.svd(realigned)
    first = (math.sqrt(singular[0]) * v[:, 0]).reshape(4, 4)
    second = (math.sqrt(singular[0]) * vh[0, :]).reshape(4, 4)
    assert np.max(np.abs(realigned - np.outer(first.ravel(), second.ravel()))) < 1e-10
    return first, second, singular


def weyl(unitary):
    decomposition = TwoQubitWeylDecomposition(unitary)
    coordinates = [float(decomposition.a), float(decomposition.b),
                   float(decomposition.c)]
    return {"angles_radians": coordinates,
            "angles_over_pi": [x / math.pi for x in coordinates],
            "nearest_pi_over_16": [round(16 * x / math.pi) for x in coordinates],
            "nearest_pi_over_16_error": [
                float(x - round(16 * x / math.pi) * math.pi / 16)
                for x in coordinates]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text())
    gates = candidate["gates"]
    edges = [(g["control"], g["target"]) for g in gates if g["gate"] == "cx"]
    assert edges[6:] == [(0, 2), (0, 2), (1, 3), (1, 3),
                         (0, 1), (0, 1), (2, 3), (2, 3)]
    first_layer = circuit_matrix(gates[13:65])
    second_layer = circuit_matrix(gates[65:])
    a, b, singular_a = factor(first_layer, (0, 2))
    c, d, singular_b = factor(second_layer, (0, 1))
    result = {"candidate": str(args.candidate),
              "layer_factorization_singular_values": {
                  "02|13": [float(x) for x in singular_a[:4]],
                  "01|23": [float(x) for x in singular_b[:4]]},
              "weyl": {"pair_02": weyl(a), "pair_13": weyl(b),
                       "pair_01": weyl(c), "pair_23": weyl(d)}}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
