"""Operator-Schmidt lower bounds for the global-parity-gauged boundary block.

Boundary chronology: RZZ(1,2,delta), inverse open-prefix CX0,2 and CX3,1,
then the first two Fourier butterflies baseline[18:38]. For a one-wire cut,
each crossing CNOT multiplies operator-Schmidt rank by at most two. If all
four one-wire cuts have rank four, every wire needs at least two CNOT
incidences, hence any exact synthesis of this block needs at least four CX.
"""

import json
import math
from pathlib import Path

import numpy as np

from check_circuit import (cnot, embedded_one_qubit, one_qubit_matrix)


ROOT = Path(__file__).resolve().parent
SPLITS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))


def matrix(delta):
    baseline = json.loads((ROOT / "baseline_18.json").read_text())["gates"]
    gates = [
        {"gate": "cx", "control": 1, "target": 2},
        {"gate": "rz", "qubit": 2, "theta": delta},
        {"gate": "cx", "control": 1, "target": 2},
        {"gate": "cx", "control": 0, "target": 2},
        {"gate": "cx", "control": 3, "target": 1},
        *baseline[18:38],
    ]
    u = np.eye(16, dtype=complex)
    for gate in gates:
        m = (cnot(gate["control"], gate["target"], 4) if gate["gate"] == "cx"
             else embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        u = m @ u
    return u


def bits(x, wires):
    result = 0
    for q in wires:
        result = (result << 1) | ((x >> (3 - q)) & 1)
    return result


def schmidt_singular_values(u, left):
    right = tuple(q for q in range(4) if q not in left)
    nleft, nright = 1 << len(left), 1 << len(right)
    realignment = np.zeros((nleft * nleft, nright * nright), dtype=complex)
    for output in range(16):
        oa, ob = bits(output, left), bits(output, right)
        for inp in range(16):
            ia, ib = bits(inp, left), bits(inp, right)
            realignment[oa * nleft + ia, ob * nright + ib] = u[output, inp]
    return np.linalg.svd(realignment, compute_uv=False)


def main():
    rows = []
    for k in range(17):
        values = {}
        u = matrix(k * math.pi / 16)
        for split in SPLITS:
            singular = schmidt_singular_values(u, split)
            rank = int(np.sum(singular > 1e-9))
            values["".join(map(str, split))] = {
                "rank": rank,
                "smallest_nonzero_singular": float(singular[rank - 1]),
            }
        rows.append({"k": k, "delta_over_pi": f"{k}/16", "splits": values,
                     "single_wire_degree_bound":
                     int(math.ceil(sum(math.ceil(math.log2(values[str(q)]["rank"]))
                                       for q in range(4)) / 2))})
    print(json.dumps({"definition": "delta=k*pi/16; cuts use physical MSB-first wires",
                      "results": rows}, indent=2))


if __name__ == "__main__":
    main()
