"""Warm-start prefix edge reroutes from the best 13-CNOT prefix deletion.

Delete exact14 prefix CNOT 4 and change one of the neighboring retained
prefix CNOT edges. The run-230 approximate local gates initialize every fit;
each edge change is a distinct 13-CNOT topology. Direct cycle loss allows
the output eigenbasis to vary. All negative results are numerical only.
"""

import argparse
import json
from pathlib import Path

from search13_multi_pair_block_rewire import fit
from search13_prefix_perturb import BLOCKS, SOURCE, decode_layers


ROOT = Path(__file__).resolve().parent
REWIRES = (
    (3, 0, 1), (3, 1, 0), (3, 1, 2), (3, 2, 1),
    (3, 2, 3), (3, 3, 2),
    (5, 0, 2), (5, 2, 0), (5, 0, 3), (5, 3, 0),
    (5, 1, 2), (5, 2, 1),
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=250)
    parser.add_argument("--seed", type=int, default=134500)
    parser.add_argument("--prefix", default="search13_prefix_rewire_warm")
    args = parser.parse_args()
    angles = decode_layers(json.loads(SOURCE.read_text()))
    rows = []
    for case, rewire in enumerate(REWIRES):
        candidate, record = fit(4, BLOCKS, args.seed + case, 0.03,
                                args.maxiter, initial_angles=angles,
                                prefix_rewire=rewire)
        path = ROOT / f"{args.prefix}_case{case}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        record.update(case=case, candidate=path.name)
        rows.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
    result = {"source": SOURCE.name, "deleted_prefix_cx": 4,
              "rewires": REWIRES, "seed": args.seed,
              "initial_perturb_sigma": 0.03, "maxiter": args.maxiter,
              "rows": rows, "best": min(rows, key=lambda row: row["loss"])}
    path = ROOT / f"{args.prefix}_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": path.name,
                      "best": result["best"]["candidate"]}), flush=True)


if __name__ == "__main__":
    main()
