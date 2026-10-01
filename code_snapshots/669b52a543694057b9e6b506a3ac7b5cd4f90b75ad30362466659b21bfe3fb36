"""Screen rank-compatible seven-CNOT block graphs with label-free fitting.

The exact 14-CNOT circuit has a fixed six-CNOT prefix and four two-CNOT
blocks. We delete one block CNOT and change at least two whole block edges.
Each retained graph has enough crossings to realize the ungauged full-tail
operator-Schmidt ranks (4 on single-wire cuts, 16 on balanced cuts). This
is only a necessary graph screen; the optimization is numerical evidence.
"""

import argparse
import json
from itertools import product
from pathlib import Path

from search13_multi_pair_block_rewire import fit


ROOT = Path(__file__).resolve().parent
PAIRS = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
BASE = ((0, 2), (1, 3), (0, 1), (2, 3))
CUTS = ((0, 1), (0, 2), (0, 3))


def graph_profile(blocks, deleted):
    edges = [edge for block in blocks for edge in (block, block)]
    edges.pop(deleted - 6)
    degrees = tuple(sum(wire in edge for edge in edges) for wire in range(4))
    crossing = tuple(sum((a in edge) != (b in edge) for edge in edges)
                     for a, b in CUTS)
    return degrees, crossing


def choose_cases(limit):
    candidates = []
    for deleted in (6, 8, 10, 12):
        for blocks in product(PAIRS, repeat=4):
            changed = sum(left != right for left, right in zip(blocks, BASE))
            if changed < 2:
                continue
            degrees, crossing = graph_profile(blocks, deleted)
            if min(degrees) < 2 or min(crossing) < 4:
                continue
            # Prefer near-baseline, degree-balanced graphs; retain diversity
            # across the four deleted block positions and the cut profile.
            score = (changed, max(degrees) - min(degrees),
                     sum(abs(c - 4) for c in crossing), blocks)
            candidates.append((score, deleted, blocks, degrees, crossing))
    chosen = []
    used_profiles = set()
    for deleted in (6, 8, 10, 12):
        pool = sorted(row for row in candidates if row[1] == deleted)
        quota = (limit + 3) // 4
        for row in pool:
            signature = (deleted, tuple(sorted(row[3])), row[4])
            if signature in used_profiles:
                continue
            chosen.append(row)
            used_profiles.add(signature)
            if sum(item[1] == deleted for item in chosen) >= quota:
                break
    return chosen[:limit], len(candidates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=12)
    parser.add_argument("--maxiter", type=int, default=120)
    parser.add_argument("--seed", type=int, default=133700)
    parser.add_argument("--prefix", default="search13_block_graph_screen")
    args = parser.parse_args()
    cases, eligible = choose_cases(args.limit)
    rows = []
    for i, (_, deleted, blocks, degrees, crossing) in enumerate(cases):
        candidate, record = fit(deleted, blocks, args.seed + i,
                                0.12, args.maxiter)
        path = ROOT / f"{args.prefix}_case{i}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        record.update(case=i, degrees=degrees, crossing=crossing,
                      candidate=path.name)
        rows.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
    result = {"base": "topology14_exact_matchgate_rational.json",
              "method": "fixed six-CNOT prefix, changed four-block graph, "
                        "direct label-free diagonalization fit",
              "eligible_graphs": eligible, "limit": args.limit,
              "maxiter": args.maxiter, "seed": args.seed,
              "rows": rows, "best": min(rows, key=lambda row: row["loss"])}
    result_path = ROOT / f"{args.prefix}_result.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": result_path.name,
                      "best": result["best"]["case"]}), flush=True)


if __name__ == "__main__":
    main()
