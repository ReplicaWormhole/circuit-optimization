"""Exact Schmidt-rank topology bounds for one shared orbit Fourier block.

Construct a label-preserving 4-qubit direct sum: Bell SWAP basis on label00,
common two-bit QFT on labels01/10, and QFT composed with Gray CX on label11.
Verify it diagonalizes the exact conjugated V4 action, compute seven exact
operator-Schmidt ranks, and enumerate the minimum edge-count multigraphs
passing the cut constraints. This is a lower bound, not CNOT synthesis.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import sympy as sp

from topology14_orbit_branch_basis import (CP_PLUS, CX01, H0, H1, IDENTITY,
                                           permutation)


ROOT = Path(__file__).resolve().parent
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
EDGES = tuple(itertools.combinations(range(4), 2))


def label(x):
    return x & 3


def position(x):
    return x >> 2


def full_block():
    fourier = H1 * CP_PLUS * H0
    bell = H0 * CX01
    branches = (bell, fourier, fourier, fourier * CX01)
    u = sp.zeros(16)
    for output in range(16):
        for inp in range(16):
            if label(output) == label(inp):
                u[output, inp] = branches[label(inp)][position(output), position(inp)]
    return u


def bit(x, q):
    return (x >> (3 - q)) & 1


def bits(x, wires):
    out = 0
    for q in wires:
        out = (out << 1) | bit(x, q)
    return out


def realign(u, left):
    right = tuple(q for q in range(4) if q not in left)
    nleft, nright = 1 << len(left), 1 << len(right)
    out = sp.zeros(nleft * nleft, nright * nright)
    for output in range(16):
        for inp in range(16):
            row = bits(output, left) * nleft + bits(inp, left)
            col = bits(output, right) * nright + bits(inp, right)
            out[row, col] = u[output, inp]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    action = json.loads((ROOT / "topology14_orbit_action_result.json").read_text())
    permutation_map = action["conjugated_permutation"]
    t = sp.zeros(16)
    for inp, output in enumerate(permutation_map):
        t[output, inp] = 1
    u = full_block()
    conjugated = sp.simplify(u * t * u.conjugate().T)
    assert all(conjugated[i, j] == 0 for i in range(16)
               for j in range(16) if i != j)
    ranks = {}
    for cut in CUTS:
        ranks["".join(map(str, cut))] = int(realign(u, cut).rank())
    bounds = [(cut, math.ceil(math.log2(ranks["".join(map(str, cut))])))
              for cut in CUTS]
    first_passing = None
    passing_count = 0
    for count in range(9):
        passing = [graph for graph in itertools.combinations_with_replacement(EDGES, count)
                   if all(sum((a in cut) != (b in cut) for a, b in graph) >= bound
                          for cut, bound in bounds)]
        if passing:
            first_passing = count
            passing_count = len(passing)
            break
    result = {"exact_diagonalization": True,
              "diagonal_labels_by_output_index": [str(conjugated[j, j])
                                                   for j in range(16)],
              "cut_ranks": ranks, "cut_based_cnot_lower_bound": first_passing,
              "passing_graph_count_at_lower_bound": passing_count,
              "scope": "one exact label-preserving branch basis; other eigenspace choices may differ"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
