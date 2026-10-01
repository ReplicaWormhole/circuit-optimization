"""Direct 13-CNOT fit after deleting one of exact14's six prefix CNOTs.

Unlike prior tail-deletion screens, this keeps the eight-CNOT matchgate
interaction graph and refits arbitrary single-qubit rotations across the
whole circuit against the label-free cycle objective. A failed fit is not
a no-go theorem.
"""

import argparse
import json
from pathlib import Path

from search13_multi_pair_block_rewire import fit


ROOT = Path(__file__).resolve().parent
BLOCKS = ((0, 2), (1, 3), (0, 1), (2, 3))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=250)
    parser.add_argument("--seed", type=int, default=134100)
    parser.add_argument("--prefix", default="search13_prefix_delete_fit")
    args = parser.parse_args()
    rows = []
    for deleted in range(6):
        for start, sigma in enumerate((0.03, 0.2)):
            seed = args.seed + 2 * deleted + start
            candidate, record = fit(deleted, BLOCKS, seed, sigma,
                                    args.maxiter)
            path = ROOT / f"{args.prefix}_p{deleted}_s{start}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            record.update(start=start, candidate=path.name)
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    result = {"base": "topology14_exact_matchgate_rational.json",
              "deleted_prefix_positions": list(range(6)),
              "tail_blocks": BLOCKS, "seed": args.seed,
              "maxiter_per_start": args.maxiter,
              "rows": rows, "best": min(rows, key=lambda row: row["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": path.name,
                      "best": result["best"]["candidate"]}), flush=True)


if __name__ == "__main__":
    main()
