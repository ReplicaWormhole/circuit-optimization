"""Exhaustively screen the eight four-mask parity models at five CNOTs.

The earlier first-pair screen covered 54 five-mask phase models. These eight
four-mask models are a distinct pi/8 phase-polynomial family. We enumerate
12**5 directed CNOT schedules and record any endpoint one CNOT from the
exact 14-CNOT circuit's six-CNOT parity prefix on the first matchgate pairs.
"""

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BASIS = (8, 4, 2, 1)
EXACT6_ROWS = (8, 5, 10, 1)


def main():
    models = json.loads((ROOT / "algebraic14_fourmask_models.json").read_text())
    supports = {frozenset(row["support"]): index for index, row in enumerate(models)}
    assert len(supports) == len(models) == 8
    targets = {}
    for a, b in ((0, 2), (2, 0), (1, 3), (3, 1)):
        rows = list(EXACT6_ROWS)
        rows[b] ^= rows[a]
        targets[tuple(rows)] = f"{a}->{b}"

    hits = []
    counts = Counter()

    def walk(rows, masks, path):
        if len(path) == 5:
            target = targets.get(rows)
            if target:
                counts[target] += 1
                for support, index in supports.items():
                    if support <= masks:
                        hits.append({"orientation": target,
                                     "model_index": index,
                                     "support": sorted(support),
                                     "path": [list(edge) for edge in path],
                                     "visited_masks": sorted(masks),
                                     "endpoint": list(rows)})
            return
        for a in range(4):
            for b in range(4):
                if a == b:
                    continue
                next_rows = list(rows)
                next_rows[b] ^= rows[a]
                mask = next_rows[b]
                next_masks = masks | {mask} if mask.bit_count() > 1 else masks
                walk(tuple(next_rows), next_masks, path + ((a, b),))

    walk(BASIS, frozenset(), ())
    result = {"phase_model_count": len(models), "path_count": 12**5,
              "all_endpoint_path_counts": dict(counts),
              "eligible_path_count": len(hits),
              "eligible_by_orientation": dict(Counter(hit["orientation"] for hit in hits)),
              "hits": hits}
    out = ROOT / "search13_firstpair_fourmask_screen_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "hits"},
                     indent=2))


if __name__ == "__main__":
    main()
