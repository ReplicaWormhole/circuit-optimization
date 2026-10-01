"""Exact modular cut-rank obstruction at a shifted exact-14 tail boundary.

Fix the six-CNOT parity prefix and the complete first F(0,2) matchgate.
The remaining suffix has six CNOTs. For each single output CNOT that joins
the terminal pairs (01)|(23), specialize its cyclotomic entries in F_193.
Nonzero modular minors certify lower bounds on exact Schmidt ranks, hence
on the CNOT count of a resynthesis of that *fixed gauged suffix*.
"""

import itertools
import json
from pathlib import Path

import numpy as np

from delete14_integer_audit import matrix


P = 193
ROOT = Path(__file__).resolve().parent
CUTS = ((0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3))
EDGES = tuple(itertools.combinations(range(4), 2))


def primitive_96_root():
    for z in range(2, P):
        if pow(z, 96, P) == 1 and all(pow(z, 96 // q, P) != 1 for q in (2, 3)):
            return z
    raise AssertionError("no primitive 96th root")


Z = primitive_96_root()
POWERS = np.array([pow(Z, k, P) for k in range(32)], dtype=np.int64)


def scalar(poly):
    return int(np.dot(POWERS, np.asarray(poly, dtype=np.int64)) % P)


def gate_matrix(gate):
    if gate["gate"] == "cx":
        out = np.zeros((16, 16), dtype=np.int64)
        c, t = gate["control"], gate["target"]
        for col in range(16):
            row = col ^ ((8 >> t) if col & (8 >> c) else 0)
            out[row, col] = 1
        return out
    q = gate["qubit"]
    local, _ = matrix(gate)
    out = np.zeros((16, 16), dtype=np.int64)
    for col in range(16):
        b = int(bool(col & (8 >> q)))
        for a in (0, 1):
            row = (col & ~(8 >> q)) | (a * (8 >> q))
            out[row, col] = scalar(local[a][b])
    return out


def realign(u, cut):
    other = tuple(q for q in range(4) if q not in cut)
    shape = [2] * 8
    order = list(cut) + [q + 4 for q in cut] + list(other) + [q + 4 for q in other]
    return u.reshape(shape).transpose(order).reshape(4 ** len(cut), 4 ** len(other))


def rank_mod(a):
    a = a.copy() % P
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i, c]), None)
        if pivot is None:
            continue
        a[[r, pivot]] = a[[pivot, r]]
        a[r] = a[r] * pow(int(a[r, c]), -1, P) % P
        for i in range(r + 1, rows):
            if a[i, c]:
                a[i] = (a[i] - a[i, c] * a[r]) % P
        r += 1
        if r == rows:
            break
    return r


def minimum_edges(rank_bits):
    crossing = np.array([[int((a in cut) != (b in cut)) for a, b in EDGES]
                         for cut in CUTS], dtype=np.int64)
    for n in range(9):
        for choice in itertools.combinations_with_replacement(range(6), n):
            counts = np.bincount(choice, minlength=6)
            if np.all(crossing @ counts >= rank_bits):
                return n, counts.tolist()
    raise AssertionError("no feasible graph within eight CNOTs")


def main():
    gates = json.loads((ROOT / "topology14_exact_matchgate_rational.json").read_text())["gates"]
    suffix = gates[31:]
    assert [(g["control"], g["target"]) for g in suffix if g["gate"] == "cx"] == [
        (1, 3), (1, 3), (0, 1), (0, 1), (2, 3), (2, 3)]
    u = np.eye(16, dtype=np.int64)
    for gate in suffix:
        u = gate_matrix(gate) @ u % P
    rows = []
    for c in range(4):
        for t in range(4):
            if c == t or (c < 2) == (t < 2):
                continue
            gauged = gate_matrix({"gate": "cx", "control": c, "target": t}) @ u % P
            ranks = [rank_mod(realign(gauged, cut)) for cut in CUTS]
            bits = [(r - 1).bit_length() for r in ranks]
            lower, witness = minimum_edges(bits)
            rows.append({"output_cnot": [c, t], "modular_ranks": ranks,
                         "cut_cnot_lower_bounds": bits,
                         "minimum_graph_edges": lower,
                         "feasible_edge_multiplicities": witness})
    result = {"field": "F_193", "primitive_96_root": Z,
              "suffix_gate_start": 31, "cuts": [list(c) for c in CUTS],
              "edge_order": [list(e) for e in EDGES], "rows": rows,
              "all_five_cnot_resyntheses_excluded": all(
                  row["minimum_graph_edges"] >= 6 for row in rows)}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
