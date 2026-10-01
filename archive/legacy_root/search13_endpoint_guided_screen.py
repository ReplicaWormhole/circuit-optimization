"""Rank five-CNOT paths one first-pair CNOT from the exact prefix endpoint.

The catalogued complete parity phase models have no eligible path at these
endpoints.  Here we seek paths visiting many of the exact14 nonlocal phase
masks, leaving an omitted mask to be absorbed jointly into the first pair.
"""

from collections import Counter
import json

from search13_joint_gauge_fit import ROOT


BASIS = (8, 4, 2, 1)
EXACT = (8, 5, 10, 1)
OLD = ((3, 1), (1, 2), (3, 2), (0, 2), (1, 2), (3, 2))
MASKS = frozenset((7, 6, 14, 11))
ORIENTATIONS = ((2, 0), (1, 3))  # first blocks with exact two-CNOT Cartan face
TARGETS = {tuple(EXACT[j] ^ EXACT[a] if j == b else EXACT[j]
                 for j in range(4)): f"{a}->{b}"
           for a, b in ORIENTATIONS}


def distance_from_deletion(path):
    return min(sum(p != q for p, q in zip(path, OLD[:deleted]+OLD[deleted+1:]))
               for deleted in range(6))


def enumerate_paths():
    rows = []
    counts = Counter()

    def walk(endpoint, masks, path):
        if len(path) == 5:
            orientation = TARGETS.get(endpoint)
            if orientation is None:
                return
            counts[orientation] += 1
            visited = set(masks)
            score = len(MASKS & visited)
            record = {"orientation": orientation, "path": path,
                      "endpoint": endpoint, "visited_masks": sorted(visited),
                      "covered_required_masks": sorted(MASKS & visited),
                      "omitted_required_masks": sorted(MASKS - visited),
                      "required_mask_count": score,
                      "deletion_distance": distance_from_deletion(path)}
            rows.append(record)
            return
        for a in range(4):
            for b in range(4):
                if a == b:
                    continue
                next_rows = list(endpoint)
                next_rows[b] ^= endpoint[a]
                new_masks = masks | {next_rows[b]} if next_rows[b].bit_count() > 1 else masks
                walk(tuple(next_rows), new_masks, path+((a, b),))

    walk(BASIS, frozenset(), ())
    return rows, counts


def main():
    rows, counts = enumerate_paths()
    eligible = [r for r in rows if r["deletion_distance"] >= 2]
    eligible.sort(key=lambda r: (-r["required_mask_count"],
                                 -r["deletion_distance"],
                                 len(r["omitted_required_masks"]), r["path"]))
    selected = {}
    for orientation in ("2->0", "1->3"):
        selected[orientation] = [r for r in eligible if r["orientation"] == orientation][:8]
    result = {"exact_endpoint": EXACT, "targets": {v: k for k, v in TARGETS.items()},
              "required_masks": sorted(MASKS), "total_paths": 12**5,
              "target_counts": dict(counts),
              "structurally_distinct_count": len(eligible),
              "coverage_histogram": dict(Counter(r["required_mask_count"] for r in eligible)),
              "selected": selected,
              "all_structurally_distinct": eligible}
    path = ROOT / "search13_endpoint_guided_screen_result.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("all_structurally_distinct", "selected")},
                     sort_keys=True))
    print(json.dumps({"selected": {k: v[:2] for k, v in selected.items()}},
                     sort_keys=True))


if __name__ == "__main__":
    main()
