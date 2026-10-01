"""Bounded exact rank screen for output row permutations within shift sectors.

This explores a gauge family distinct from one/two row swaps or Hadamard
blocks. Modular rank is a rigorous lower bound on characteristic-zero rank;
sampling cannot rule out the unsampled permutations.
"""

import argparse
from collections import Counter
from itertools import permutations
import json
from pathlib import Path

import numpy as np

from search13_fulltail_invariant_gauge_rank import setup
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment
from search13_sparse_gauge_clifford_pairs import target_mod


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_sector_row_permutations_result.json"
MASKS = (83, 163)
PRIME = 97


def rank_pair(base, perm):
    target = base[perm, :]
    return tuple(rank_mod(realignment(target, m), PRIME) for m in MASKS)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=30113)
    parser.add_argument("--samples", type=int, default=30000)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    base, _, _ = target_mod(PRIME)
    _, _, _, sectors = setup()
    rng = np.random.default_rng(args.seed)
    rows = np.arange(16)
    best = (99, None)
    profiles = Counter()
    for index in range(args.samples):
        perm = rows.copy()
        for sector in sectors:
            perm[np.asarray(sector)] = rng.permutation(sector)
        ranks = rank_pair(base, perm)
        profiles[ranks] += 1
        score = sum(ranks)
        if score < best[0]:
            best = (score, {"sample_index": index, "row_permutation": perm.tolist(),
                            "mod97_ranks": ranks})
            print(json.dumps({"sample": index, "score": score, "ranks": ranks}), flush=True)
    result = {"source": "exact run266 gauge times exact14 tail",
              "family": "independent output row permutations within each shift eigenspace",
              "sector_sizes": list(map(len, sectors)), "seed": args.seed,
              "samples": args.samples, "mask_order": MASKS,
              "mod97_rank_profiles": {f"{a},{b}": n for (a, b), n in sorted(profiles.items())},
              "best": best[1],
              "scope": "random sample; modular ranks lower bound exact ranks; not exhaustive"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"best": best[1], "profiles": len(profiles)}), flush=True)


if __name__ == "__main__":
    main()
