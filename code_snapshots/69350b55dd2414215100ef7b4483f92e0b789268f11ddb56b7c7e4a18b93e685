"""Rank-screen 14-CNOT topology proposals for a distinct exact eigenbasis.

The input target is the parity-conditioned fermionic Fourier eigenbasis.
Operator-Schmidt cut bounds are necessary for synthesis of that *specific*
unitary, and may discard a topology before numerical local-gate fitting.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CUTS = [(0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3)]


def schmidt_rank(matrix, cut):
    other = [q for q in range(4) if q not in cut]
    axes = list(cut) + [q + 4 for q in cut] + other + [q + 4 for q in other]
    reshaped = matrix.reshape([2] * 8).transpose(axes).reshape(
        4 ** len(cut), 4 ** len(other))
    singular = np.linalg.svd(reshaped, compute_uv=False)
    return int(np.count_nonzero(singular > 1e-9)), float(
        singular[singular > 1e-9][-1])


def graph_passes(edges, bounds):
    return all(sum((a in cut) != (b in cut) for a, b in edges) >= bound
               for cut, bound in zip(CUTS, bounds))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    target = np.load(args.target)
    if target.shape != (16, 16):
        raise ValueError("expected a 16x16 target")
    ranks = [schmidt_rank(target, cut) for cut in CUTS]
    bounds = [int(math.ceil(math.log2(rank))) for rank, _ in ranks]
    topologies = {}
    for name in ("delete14_raw_top.json", "delete14_topology2_raw_top.json"):
        circuit = json.loads((ROOT / name).read_text())
        edges = tuple((gate["control"], gate["target"])
                      for gate in circuit["gates"] if gate["gate"] == "cx")
        topologies.setdefault(edges, []).append({"source": name, "loss": None})
    for name in ("delete14_anneal_exact_result.json",
                 "delete14_anneal_fivemask2_result.json"):
        batch = json.loads((ROOT / name).read_text())
        for row in batch["records"]:
            edges = tuple(map(tuple, row["topology"]))
            topologies.setdefault(edges, []).append(
                {"source": name, "proposal": row["proposal"],
                 "loss": row["optimized_loss"]})
    records = []
    for edges, sources in topologies.items():
        record = {"edges": [list(edge) for edge in edges], "sources": sources,
                  "rank_compatible": graph_passes(edges, bounds)}
        records.append(record)
    result = {"target": str(args.target),
              "cut_order": [list(cut) for cut in CUTS],
              "ranks": [rank for rank, _ in ranks],
              "smallest_nonzero_singular": [value for _, value in ranks],
              "cnot_crossing_lower_bounds": bounds,
              "topologies_screened": len(records),
              "rank_compatible": sum(row["rank_compatible"] for row in records),
              "records": records}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
