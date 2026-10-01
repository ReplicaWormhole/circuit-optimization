"""Bounded two-eigenspace-gauge, interleaved seven-CNOT tail search.

The exact six-CNOT prefix is retained.  Two disjoint real Givens rotations
within distinct degenerate eigenspaces act on the full eight-CNOT tail.
Arbitrary local SU(2) gates are then fitted jointly across both first-pair
and final-pair interactions.  This is a numerical screen only.
"""

import argparse
import json

import numpy as np

from check_circuit import evaluate
from delete14_search import fixed_matrix
from search13_eigenrow_gauge_seven_fit import BASE, ROOT, fit


GAUGE_PAIRS = ((8, 11), (1, 14))  # eigenvalues +1 and -i respectively
SCHEDULES = (
    ((0, 2), (0, 1), (0, 2), (0, 3), (2, 3), (1, 3), (0, 1)),
    ((0, 2), (0, 1), (1, 3), (0, 3), (0, 2), (2, 3), (0, 1)),
)


def target_and_prefix():
    gates = json.loads(BASE.read_text())["gates"]
    prefix, tail = gates[:13], gates[13:]
    assert sum(g["gate"] == "cx" for g in prefix) == 6
    assert sum(g["gate"] == "cx" for g in tail) == 8
    target = np.eye(16, dtype=complex)
    for gate in tail:
        target = fixed_matrix(gate) @ target
    for first, second in GAUGE_PAIRS:
        a, b = target[first].copy(), target[second].copy()
        target[first] = (a-b)/np.sqrt(2)
        target[second] = (a+b)/np.sqrt(2)
    return prefix, target


def balanced_ranks(target):
    tensor = target.reshape((2,)*8)
    out = {}
    for cut in ((0, 1), (0, 2), (0, 3)):
        # Each cut acts on both output and input wire indices.
        left = tuple(cut) + tuple(4+j for j in cut)
        right = tuple(j for j in range(8) if j not in left)
        mat = tensor.transpose(left+right).reshape(16, 16)
        out["".join(map(str, cut))] = int(np.linalg.matrix_rank(mat, tol=1e-10))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=350)
    parser.add_argument("--seed", type=int, default=25010)
    args = parser.parse_args()
    prefix, target = target_and_prefix()
    ranks = balanced_ranks(target)
    rows = []
    for index, edges in enumerate(SCHEDULES):
        capacities = {"".join(map(str, cut)): 2**sum(
            (a in cut) != (b in cut) for a, b in edges)
            for cut in ((0, 1), (0, 2), (0, 3))}
        assert all(ranks[k] <= capacities[k] for k in ranks)
        tail, row = fit(edges, args.seed+index, args.maxiter, target)
        candidate = {"n": 4, "gates": prefix+tail}
        path = ROOT / f"search13_joint_gauge_candidate_{index}.json"
        path.write_text(json.dumps(candidate, indent=2)+"\n")
        checked = evaluate(candidate)
        row.update(schedule_index=index, edges=edges, cut_capacities=capacities,
                   candidate=path.name,
                   off_diagonal_error=checked["off_diagonal_error"],
                   valid_diagonalizer=checked["valid_diagonalizer"])
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    result = {"gauge_pairs": GAUGE_PAIRS, "angle": "pi/4 real Givens",
              "balanced_ranks": ranks, "schedules": SCHEDULES,
              "maxiter": args.maxiter, "seed": args.seed, "rows": rows}
    out = ROOT / "search13_joint_gauge_fit_result.json"
    out.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
