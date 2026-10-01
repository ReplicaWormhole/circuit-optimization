"""Exact orbit-wise Fourier block and cut ranks after one shear.

Assign output Fourier mode k to the k-th state of each transformed V4 orbit.
This is one canonical exact eigenbasis, not an optimization over the large
degenerate-eigenspace freedom. Compute all seven operator-Schmidt ranks and
the cut-based CNOT lower bound for this chosen block.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import sympy as sp

from topology14_orbit_block_rank import CUTS, EDGES, realign


ROOT = Path(__file__).resolve().parent


def cycles(permutation):
    unseen = set(range(16))
    result = []
    while unseen:
        start = min(unseen)
        orbit = [start]
        current = permutation[start]
        while current != start:
            orbit.append(current)
            current = permutation[current]
        unseen.difference_update(orbit)
        result.append(orbit)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    action = json.loads((ROOT / "topology14_one_shear_action_result.json").read_text())
    permutation = action["induced_permutation"]
    orbit_list = cycles(permutation)
    u = sp.zeros(16)
    for orbit in orbit_list:
        length = len(orbit)
        root = sp.I if length == 4 else (-1 if length == 2 else 1)
        norm = sp.sqrt(length)
        for k, output in enumerate(orbit):
            for j, inp in enumerate(orbit):
                u[output, inp] = root ** (k * j) / norm
    t = sp.zeros(16)
    for inp, output in enumerate(permutation):
        t[output, inp] = 1
    conjugated = sp.simplify(u * t * u.conjugate().T)
    assert all(conjugated[row, col] == 0 for row in range(16)
               for col in range(16) if row != col)
    ranks = {"".join(map(str, cut)): int(realign(u, cut).rank()) for cut in CUTS}
    bounds = [(cut, math.ceil(math.log2(ranks["".join(map(str, cut))])))
              for cut in CUTS]
    minimum = None
    for count in range(11):
        if any(all(sum((a in cut) != (b in cut) for a, b in graph) >= bound
                   for cut, bound in bounds)
               for graph in itertools.combinations_with_replacement(EDGES, count)):
            minimum = count
            break
    result = {"orbits": orbit_list,
              "exact_diagonalization": True,
              "output_eigenvalues": [str(conjugated[j, j]) for j in range(16)],
              "cut_ranks": ranks,
              "cut_based_cnot_lower_bound": minimum,
              "scope": "canonical orbit-wise Fourier basis; no eigenbasis optimization"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
