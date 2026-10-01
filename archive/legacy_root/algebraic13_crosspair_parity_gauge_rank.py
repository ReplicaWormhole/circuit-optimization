"""Exact modular cut-rank screen of cross-pair parity gauges on F01*F23.

For every parity mask meeting both terminal pairs and k=0,...,15, form
G=exp(i*k*pi/16*Z_mask) and S=F01 tensor F23, where
F=exp(i*pi*(XX+YY)/8). A nonzero rank after specialization in a prime
field proves at least that rank over the cyclotomic field. The graph test
gives necessary CNOT counts for synthesis of this fixed gauged suffix.
Angles k+16 differ only by a global minus sign. Products of multiple parity
gauges and a shifted circuit boundary are outside this finite scan.
"""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PRIMES = (193, 257)
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
EDGES = tuple(itertools.combinations(range(4), 2))
MASKS = tuple(m for m in range(1, 16) if m & 0b1100 and m & 0b0011)


def primitive_32_root(p):
    for z in range(2, p):
        if pow(z, 32, p) == 1 and pow(z, 16, p) != 1:
            return z
    raise AssertionError(f"no primitive 32nd root in F_{p}")


def suffix_mod(p, z):
    z8 = pow(z, 4, p)
    rt2 = (z8 + pow(z8, -1, p)) % p
    assert rt2 * rt2 % p == 2
    a = pow(rt2, -1, p)
    ia = pow(z, 8, p) * a % p
    f = np.array([[1, 0, 0, 0], [0, a, ia, 0],
                  [0, ia, a, 0], [0, 0, 0, 1]], dtype=np.int64)
    return np.kron(f, f) % p


def gauged_suffix(p, z, suffix, mask, k):
    phases = [pow(z, k if (basis & mask).bit_count() % 2 == 0 else -k, p)
              for basis in range(16)]
    return (np.asarray(phases, dtype=np.int64)[:, None] * suffix) % p


def realign(unitary, cut):
    other = tuple(q for q in range(4) if q not in cut)
    order = list(cut) + [q + 4 for q in cut] + list(other) + [q + 4 for q in other]
    return unitary.reshape([2] * 8).transpose(order).reshape(
        4 ** len(cut), 4 ** len(other))


def rank_mod(matrix, p):
    a = matrix.copy() % p
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r, col]), None)
        if pivot is None:
            continue
        a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, p) % p
        for row in range(rank + 1, rows):
            if a[row, col]:
                a[row] = (a[row] - a[row, col] * a[rank]) % p
        rank += 1
        if rank == rows:
            break
    return rank


def minimum_graph_edges(rank_bits):
    crossing = np.array([[int((a in cut) != (b in cut)) for a, b in EDGES]
                         for cut in CUTS], dtype=np.int64)
    for n in range(9):
        for choice in itertools.combinations_with_replacement(range(6), n):
            counts = np.bincount(choice, minlength=6)
            if np.all(crossing @ counts >= rank_bits):
                return n, counts.tolist()
    raise AssertionError("no cut-compatible graph within eight CNOTs")


def main():
    fields = []
    per_field = {}
    for p in PRIMES:
        z = primitive_32_root(p)
        s = suffix_mod(p, z)
        fields.append({"prime": p, "primitive_32_root": z})
        for mask in MASKS:
            for k in range(16):
                u = gauged_suffix(p, z, s, mask, k)
                per_field[(p, mask, k)] = [rank_mod(realign(u, cut), p)
                                           for cut in CUTS]
    rows = []
    for mask in MASKS:
        for k in range(16):
            ranks_by_field = [per_field[(p, mask, k)] for p in PRIMES]
            certified = [max(ranks) for ranks in zip(*ranks_by_field)]
            bits = [(rank - 1).bit_length() for rank in certified]
            graph_min, witness = minimum_graph_edges(np.array(bits))
            rows.append({"mask": mask, "wires": [q for q in range(4)
                                                   if mask & (8 >> q)],
                         "angle_numerator_pi_over_16": k,
                         "modular_ranks_by_field": ranks_by_field,
                         "certified_rank_lower_bounds": certified,
                         "cut_cnot_lower_bounds": bits,
                         "minimum_graph_edges": graph_min,
                         "feasible_edge_multiplicities": witness})
    result = {"suffix": "F01 tensor F23, F=exp(i*pi*(XX+YY)/8)",
              "gauge": "exp(i*k*pi/16*Z_mask)",
              "fields": fields, "cuts": [list(c) for c in CUTS],
              "edge_order": [list(e) for e in EDGES],
              "sampled_masks": list(MASKS), "sampled_k": list(range(16)),
              "rows": rows,
              "three_cnot_compatible_cases": [
                  {"mask": r["mask"], "k": r["angle_numerator_pi_over_16"]}
                  for r in rows if r["minimum_graph_edges"] <= 3]}
    path = ROOT / "algebraic13_crosspair_parity_gauge_rank_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "fields": fields,
                      "cases": len(rows),
                      "three_cnot_compatible_cases": result[
                          "three_cnot_compatible_cases"],
                      "min_graph_edges_histogram": {str(n): sum(
                          r["minimum_graph_edges"] == n for r in rows)
                          for n in range(9)}}, indent=2))


if __name__ == "__main__":
    main()
