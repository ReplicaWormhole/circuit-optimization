"""Enumerate short output relabelings exposing a local degenerate eigenspace.

For every output CNOT sequence of bounded length, test all wire pairs.
The exchange type mixes |01> with |10>; the even type mixes |00> with |11>.
This is a combinatorial screen only, not a circuit synthesis or CNOT bound.
"""

import argparse
import itertools
import json
from pathlib import Path

from exact_check import exact_eigenvalue_labels


ROOT = Path(__file__).resolve().parent
EDGES = [(c, t) for c in range(4) for t in range(4) if c != t]


def cx_index(x, edge):
    c, t = edge
    return x ^ (1 << (3 - t)) if x & (1 << (3 - c)) else x


def relabel(labels, edges):
    return [labels[apply_inverse(y, edges)] for y in range(16)]


def apply_inverse(y, edges):
    for edge in reversed(edges):
        y = cx_index(y, edge)
    return y


def compatible(labels, pair, kind):
    q, r = pair
    bits = (1 << (3 - q), 1 << (3 - r))
    for spectator in range(16):
        if spectator & (bits[0] | bits[1]):
            continue
        a, b = ((bits[0], bits[1]) if kind == "exchange" else
                (0, bits[0] | bits[1]))
        if labels[spectator | a] != labels[spectator | b]:
            return False
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--max-relabel-cx", type=int, default=3)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if not 0 <= args.max_relabel_cx <= 3:
        raise ValueError("bounded screen supports 0 to 3 output CNOTs")
    base = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    labels = exact_eigenvalue_labels(base)["output_labels"]
    records = []
    for length in range(args.max_relabel_cx + 1):
        for edges in itertools.product(EDGES, repeat=length):
            if any(edges[i] == edges[i + 1] for i in range(length - 1)):
                continue
            shifted = relabel(labels, edges)
            for pair in itertools.combinations(range(4), 2):
                for kind in ("exchange", "even"):
                    if compatible(shifted, pair, kind):
                        records.append({"output_cx": [list(edge) for edge in edges],
                                        "wire_pair": list(pair), "kind": kind,
                                        "relabeled_labels": shifted})
    result = {"max_relabel_cx": args.max_relabel_cx,
              "compatible_count": len(records), "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    counts = {}
    for rec in records:
        key = (len(rec["output_cx"]), rec["kind"])
        counts[str(key)] = counts.get(str(key), 0) + 1
    print(json.dumps({"compatible_count": len(records), "counts": counts,
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
