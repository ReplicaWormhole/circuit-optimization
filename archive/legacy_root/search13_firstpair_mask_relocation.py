"""Exact necessary screen for moving one parity mask into the first pair block.

For each catalogued pi/8 parity model, visit all but one required nonlocal
mask with five directed CNOTs. The omitted mask must become ZZ on 02 or 13
at the five-CNOT endpoint. The endpoint mismatch to the exact six-CNOT
prefix must act on that same pair alone, so it can be included in the first
two-qubit block without changing the other pair or paying another CNOT.
"""

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BASIS = (8, 4, 2, 1)
EXACT6_ROWS = (8, 5, 10, 1)
PAIRS = ((0, 2), (1, 3))


def pair_local_endpoints():
    targets = {}
    for pair in PAIRS:
        a, b = pair
        other = {0, 1, 2, 3} - set(pair)
        span = (EXACT6_ROWS[a], EXACT6_ROWS[b], EXACT6_ROWS[a] ^ EXACT6_ROWS[b])
        for row_a in span:
            for row_b in span:
                if row_a == row_b:
                    continue
                rows = list(EXACT6_ROWS)
                rows[a], rows[b] = row_a, row_b
                rows = tuple(rows)
                assert all(rows[j] == EXACT6_ROWS[j] for j in other)
                targets.setdefault(rows, []).append(pair)
    return targets


def main():
    models = []
    for family in ("fourmask", "fivemask"):
        rows = json.loads((ROOT / f"algebraic14_{family}_models.json").read_text())
        for index, record in enumerate(rows):
            models.append((family, index, frozenset(record["support"]),
                           record["coefficients_pi_over_8"]))
    assert len(models) == 62
    targets = pair_local_endpoints()
    counts = Counter()
    hits = []

    def walk(rows, seen, path):
        if len(path) == 5:
            pairs = targets.get(rows)
            if pairs is None:
                return
            counts["pair_local_paths"] += 1
            for pair in pairs:
                if rows == EXACT6_ROWS:
                    counts["identity_endpoint_paths"] += 1
                pair_mask = rows[pair[0]] ^ rows[pair[1]]
                for family, index, support, coeffs in models:
                    missing = support - seen
                    if len(missing) != 1:
                        continue
                    counts["exactly_one_missing"] += 1
                    mask = next(iter(missing))
                    if mask != pair_mask:
                        continue
                    counts["pair_zz_relocations"] += 1
                    hits.append({"family": family, "model_index": index,
                                 "support": sorted(support),
                                 "missing_mask": mask,
                                 "missing_coefficient_pi_over_8": coeffs[str(mask)],
                                 "pair": list(pair), "endpoint": list(rows),
                                 "path": [list(edge) for edge in path],
                                 "visited_nonlocal_masks": sorted(seen)})
            return
        for a in range(4):
            for b in range(4):
                if a == b:
                    continue
                next_rows = list(rows)
                next_rows[b] ^= rows[a]
                mask = next_rows[b]
                next_seen = seen | {mask} if mask.bit_count() > 1 else seen
                walk(tuple(next_rows), next_seen, path + ((a, b),))

    walk(BASIS, frozenset(), ())
    result = {"models": len(models), "ordered_five_cnot_paths": 12**5,
              "pair_local_endpoints": len(targets),
              "counts": dict(counts), "hits": hits}
    path = ROOT / "search13_firstpair_mask_relocation_result.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "hits"},
                     indent=2))


if __name__ == "__main__":
    main()
