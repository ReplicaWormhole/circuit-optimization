"""Bounded 13-CNOT search: delete in one matchgate block, reroute an adjacent block.

The surviving CNOT in the deleted block retains its exact-14 edge.  A CNOT
in a neighboring block is rerouted to the deleted block's edge.  Each case
gets two independent local SU(2)-layer fits using the label-free loss from
search13_reroute_neighbor.py.  This is numerical evidence only.
"""

import argparse
import json
from pathlib import Path

from search13_reroute_neighbor import optimize


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "topology14_exact_matchgate_rational.json"
CASES = [
    (6, 8, (0, 2)),
    (8, 6, (1, 3)),
    (8, 10, (1, 3)),
    (10, 8, (0, 1)),
    (10, 12, (0, 1)),
    (12, 10, (2, 3)),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--seed", type=int, default=13250)
    parser.add_argument("--prefix", type=str,
                        default="search13_reroute_adjacent")
    args = parser.parse_args()
    results = []
    for case_index, (deleted, rerouted, edge) in enumerate(CASES):
        for start in range(2):
            seed = args.seed + 2 * case_index + start
            sigma = 0.04 if start == 0 else 0.3
            candidate, record = optimize(deleted, rerouted, edge, seed,
                                         sigma, args.maxiter)
            path = ROOT / f"{args.prefix}_d{deleted}_r{rerouted}_{edge[0]}{edge[1]}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            record["candidate"] = path.name
            results.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    summary = {"base": BASE.name, "results": results,
               "best": min(results, key=lambda r: r["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"result_file": path.name, "best": summary["best"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
