"""Screen omitted parity phases that become local on the first pair.

This complements search13_firstpair_mask_relocation.py: a missing mask can
move into a two-qubit block's local frame if it equals one endpoint row.
"""

import json
from collections import Counter
from pathlib import Path

from search13_firstpair_mask_relocation import BASIS, ROOT, pair_local_endpoints


def main():
    models = []
    for family in ("fourmask", "fivemask"):
        data = json.loads((ROOT / f"algebraic14_{family}_models.json").read_text())
        for index, row in enumerate(data):
            models.append((family, index, frozenset(row["support"]),
                           row["coefficients_pi_over_8"]))
    targets = pair_local_endpoints()
    hits = []
    counts = Counter()

    def walk(rows, seen, path):
        if len(path) == 5:
            pairs = targets.get(rows)
            if pairs is None:
                return
            for pair in pairs:
                for family, index, support, coeffs in models:
                    missing = support - seen
                    if len(missing) != 1:
                        continue
                    mask = next(iter(missing))
                    if mask not in (rows[pair[0]], rows[pair[1]]):
                        continue
                    counts[family] += 1
                    hits.append({"family": family, "model_index": index,
                                 "support": sorted(support),
                                 "missing_mask": mask,
                                 "missing_coefficient_pi_over_8": coeffs[str(mask)],
                                 "pair": list(pair),
                                 "local_wire": pair[0] if mask == rows[pair[0]] else pair[1],
                                 "endpoint": list(rows),
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
              "local_relocations": len(hits), "by_family": dict(counts),
              "hits": hits}
    out = ROOT / "search13_firstpair_mask_local_screen_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "hits"}, indent=2))


if __name__ == "__main__":
    main()
