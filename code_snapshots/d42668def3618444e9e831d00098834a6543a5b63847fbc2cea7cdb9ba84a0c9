"""Extract two-qubit Weyl invariants from the numerical 14-CX tail.

The eight tail CNOTs occur as four same-pair doublets. For each doublet,
only the local layer between its CNOTs controls its nonlocal Weyl class;
one-wire gates before and after are local gauges. Compute KAK coordinates
of CX (A tensor B) CX and compare to simple pi fractions.
"""

import json
import math
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator
from qiskit.synthesis import TwoQubitWeylDecomposition

from check_circuit import one_qubit_matrix
from topology14_eightcx_opt import layers_from_suffix


ROOT = Path(__file__).parent


def cx_two():
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 0, 1],
                     [0, 0, 1, 0]], dtype=complex)


def main():
    candidate = json.loads(
        (ROOT / "algebraic14_pure_orbit_refined_candidate.json").read_text())
    suffix = candidate["gates"][13:]
    angles, edges = layers_from_suffix(suffix)
    assert edges == [(0, 2), (0, 2), (1, 3), (1, 3),
                     (0, 1), (0, 1), (2, 3), (2, 3)]
    # Extract middle-layer local matrices directly from chronological gates,
    # avoiding any Euler-parameter convention issues.
    cx_seen = 0
    middle = {1: [], 3: [], 5: [], 7: []}
    for gate in suffix:
        if gate["gate"] == "cx":
            cx_seen += 1
        elif cx_seen in middle:
            middle[cx_seen].append(gate)
    rows = []
    cx = cx_two()
    for pair_index, edge in enumerate(edges[::2]):
        layer = middle[2 * pair_index + 1]
        mats = [np.eye(2, dtype=complex), np.eye(2, dtype=complex)]
        for gate in layer:
            if gate["qubit"] in edge:
                local_index = edge.index(gate["qubit"])
                mats[local_index] = one_qubit_matrix(gate) @ mats[local_index]
        # All doublets have physical control=first endpoint.
        block = cx @ np.kron(mats[0], mats[1]) @ cx
        kak = TwoQubitWeylDecomposition(block)
        coords = [float(kak.a), float(kak.b), float(kak.c)]
        rows.append({
            "pair": edge,
            "middle_layer_index": 2 * pair_index + 1,
            "weyl_coordinates": coords,
            "weyl_over_pi": [value / math.pi for value in coords],
            "nearest_pi_over_32": [round(32 * value / math.pi)
                                   for value in coords],
            "deviation_from_pi_over_32": [
                abs(value - round(32 * value / math.pi)
                    * math.pi / 32) for value in coords],
            "middle_local_matrix_0": [
                [str(value) for value in row] for row in mats[0]],
            "middle_local_matrix_1": [
                [str(value) for value in row] for row in mats[1]],
        })
    result = {"cnot_count": 14, "pair_blocks": rows}
    (ROOT / "algebraic14_pair_weyl_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: [{k: value for k, value in row.items()
                              if not k.startswith("middle_local_matrix")}
                             for row in rows] if key == "pair_blocks" else result[key]
                      for key in result}, indent=2))


if __name__ == "__main__":
    main()
