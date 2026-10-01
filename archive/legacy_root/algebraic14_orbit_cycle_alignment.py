"""Check position-cycle alignment for all two-shear orbit coordinates."""

import json
from collections import Counter
from pathlib import Path

from algebraic14_orbit_shear import orbits, parity, plane_direction, shear
from algebraic14_orbit_two_shears import form_cost, unique_shears


ROOT = Path(__file__).parent


def cycle_signature(cycle, direction):
    u, v = direction[1:3]
    anchor = min(cycle)
    local = {anchor ^ diff: index for index, diff in
             enumerate((0, u, v, u ^ v))}
    indices = [local[x] for x in cycle]
    offset = indices.index(0)
    return tuple(indices[offset:] + indices[:offset])


def main():
    cycles = [cycle for cycle in orbits() if len(cycle) == 4]
    items = list(unique_shears().items())
    histogram = Counter()
    records = []
    for first_perm, first_form in items:
        intermediate = [[first_perm[x] for x in cycle] for cycle in cycles]
        for second_perm, second_form in items:
            transformed = [[second_perm[x] for x in cycle]
                           for cycle in intermediate]
            directions = [plane_direction(cycle)
                          for cycle in transformed]
            if directions[0] is None or directions[0] != directions[1] or directions[0] != directions[2]:
                continue
            signatures = [cycle_signature(cycle, directions[0])
                          for cycle in transformed]
            multiplicity = max(Counter(signatures).values())
            histogram[multiplicity] += 1
            records.append({"first": first_form, "second": second_form,
                            "direction": directions[0],
                            "cycle_signatures": signatures,
                            "max_identical_cycle_count": multiplicity,
                            "score": form_cost(first_form)[0] +
                                     form_cost(second_form)[0],
                            "transformed_orbits": transformed})
    records.sort(key=lambda record: (-record["max_identical_cycle_count"],
                                     record["score"],
                                     record["first"], record["second"]))
    result = {"solutions": len(records),
              "max_identical_cycle_histogram": dict(histogram),
              "best": records[0] if records else None,
              "top_records": records[:16]}
    (ROOT / "algebraic14_orbit_cycle_alignment_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("solutions", "max_identical_cycle_histogram", "best")},
                     indent=2))


if __name__ == "__main__":
    main()
