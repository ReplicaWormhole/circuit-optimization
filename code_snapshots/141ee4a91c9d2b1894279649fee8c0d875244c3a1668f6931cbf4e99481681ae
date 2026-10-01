"""Bounded seven-CNOT fits on three exact-rank-compatible tail schedules.

Fit arbitrary SU(2) gates between all CNOTs first to the exact gauged
tail by process overlap, then to the label-free cycle. For the closest
schedule, also warm start the direct fit from a deleted/rewired exact8 tail.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from check_circuit import evaluate, one_qubit_matrix
from search13_eigenrow_gauge_seven_fit import fit as process_fit
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
WITNESS = ROOT / "search13_fulltail_invariant_exact_witness_result.json"
OUT = ROOT / "search13_fulltail_rank8_fit_result.json"
SCHEDULES = (
    ((0, 2), (0, 2), (0, 3), (1, 3), (0, 1), (0, 1), (2, 3)),
    ((0, 2), (0, 3), (1, 2), (1, 3), (0, 1), (0, 3), (2, 3)),
    ((0, 2), (0, 3), (1, 2), (1, 3), (0, 1), (0, 3), (1, 2)),
)


def exact_gauge_numeric():
    witness = json.loads(WITNESS.read_text())
    zeta = sp.exp(2 * sp.pi * sp.I / 96)
    g = np.zeros((16, 16), dtype=complex)
    for row, columns in witness["gauge_sparse_rows"].items():
        for col, expression in columns.items():
            g[int(row), int(col)] = complex(sp.N(sp.sympify(
                expression, locals={"zeta": zeta}), 17))
    assert np.max(np.abs(g @ g.conj().T - np.eye(16))) < 1e-12
    return g


def matrix_to_axis(m):
    a, b = m[0, 0], m[0, 1]
    components = np.array([-b.imag, -b.real, -a.imag])
    sine = np.linalg.norm(components)
    if sine < 1e-13:
        return np.zeros(3)
    angle = 2 * math.atan2(sine, a.real)
    return components * (angle / sine)


def baseline_warm_angles():
    """Delete the last old tail CNOT and rewire old tail CNOT 2 to 03."""
    source = json.loads(BASE.read_text())["gates"][13:]
    chunks = [[]]
    old_cx = -1
    for gate in source:
        if gate["gate"] == "cx":
            old_cx += 1
            if old_cx != 7:
                chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 8 and old_cx == 7
    result = np.zeros((8, 4, 3))
    for slot, chunk in enumerate(chunks):
        for qubit in range(4):
            m = np.eye(2, dtype=complex)
            for gate in chunk:
                if gate["qubit"] == qubit:
                    m = one_qubit_matrix(gate) @ m
            result[slot, qubit] = matrix_to_axis(m)
    return result


def save_candidate(prefix_gates, tail_gates, name):
    candidate = {"n": 4, "gates": prefix_gates + tail_gates}
    path = ROOT / name
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    return path.name, evaluate(candidate)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process-maxiter", type=int, default=250)
    parser.add_argument("--direct-maxiter", type=int, default=350)
    parser.add_argument("--seed", type=int, default=26700)
    args = parser.parse_args()
    tail, prefix_matrix, _, _ = setup()
    target = exact_gauge_numeric() @ tail
    prefix_gates = json.loads(BASE.read_text())["gates"][:13]
    rows = []
    for case, edges in enumerate(SCHEDULES):
        gates, process = process_fit(edges, args.seed + case,
                                     args.process_maxiter, target)
        path, check = save_candidate(prefix_gates, gates,
              f"search13_fulltail_rank8_fit_case{case}_process.json")
        process.update(candidate=path, off_diagonal_error=check["off_diagonal_error"],
                       valid_diagonalizer=check["valid_diagonalizer"])
        angles = decode_local_layers(gates)
        direct_gates, direct = direct_fit(edges, prefix_matrix, angles,
                                         args.direct_maxiter)
        path, check = save_candidate(prefix_gates, direct_gates,
              f"search13_fulltail_rank8_fit_case{case}_direct.json")
        direct.update(candidate=path, off_diagonal_error=check["off_diagonal_error"],
                      valid_diagonalizer=check["valid_diagonalizer"])
        row = {"case": case, "edges": edges, "process": process,
               "direct": direct}
        if case == 0:
            warm_gates, warm = direct_fit(edges, prefix_matrix,
                    baseline_warm_angles(), args.direct_maxiter)
            path, check = save_candidate(prefix_gates, warm_gates,
                  "search13_fulltail_rank8_fit_case0_baseline_warm.json")
            warm.update(candidate=path,
                        off_diagonal_error=check["off_diagonal_error"],
                        valid_diagonalizer=check["valid_diagonalizer"])
            row["baseline_warm_direct"] = warm
        rows.append(row)
        print(json.dumps({"case": case, "process_loss": process["loss"],
                          "direct_loss": direct["loss"],
                          "baseline_warm_loss": row.get("baseline_warm_direct", {})
                                                 .get("loss")}), flush=True)
    result = {"exact_gauge_source": WITNESS.name,
              "source": BASE.name, "schedules": SCHEDULES,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter, "seed": args.seed,
              "rows": rows, "scope": "three graph/orders, one process and direct fit each; one extra baseline warm direct fit; numerical synthesis only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
