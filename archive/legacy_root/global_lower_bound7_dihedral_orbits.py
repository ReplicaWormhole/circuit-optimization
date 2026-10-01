"""Canonical D4 classes of six-CNOT schedules surviving proven filters.

Only physical dihedral relabelings are quotiented. We do not identify
chronological reorderings, CNOT directions, or arbitrary wire permutations.
Survival is a necessary-condition status, never a circuit witness.
"""

from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path

import global_lower_bound7_schedule_filter as prior
from global_lower_bound7_adjacent_pair_support_filter import first_commuting_cone


ROOT = Path(__file__).resolve().parent


def adjacent_pair_rejection(schedule):
    for j in range(4):
        cone, _ = first_commuting_cone(schedule, j)
        if len(cone) == 2 and tuple(sorted(cone)) not in prior.OPPOSITE:
            return True
    return False


def survives(schedule):
    deg = prior.degrees(schedule)
    return (min(deg) >= 2 and prior.all_forward_cones_full(schedule)
            and prior.tail_rejection(schedule, deg) is None
            and not adjacent_pair_rejection(schedule))


def d4_images(schedule):
    for reflection in (False, True):
        for offset in range(4):
            transformed = []
            for a, b in schedule:
                mapped = (offset - a, offset - b) if reflection else (
                    offset + a, offset + b)
                transformed.append(tuple(sorted(x % 4 for x in mapped)))
            yield tuple(transformed)


def canonical(schedule):
    return min(d4_images(schedule))


def active_profiles(schedule):
    deg = prior.degrees(schedule)
    degree_two = tuple(sorted(len(prior.selected_tail(schedule, j)[2])
                              for j in range(4) if deg[j] == 2))
    first_cones = tuple(sorted(len(first_commuting_cone(schedule, j)[0])
                               for j in range(4)))
    return degree_two, first_cones


def main():
    classes = defaultdict(list)
    survivor_count = 0
    for schedule in product(prior.EDGES, repeat=6):
        if not survives(schedule):
            continue
        survivor_count += 1
        classes[canonical(schedule)].append(schedule)
    assert survivor_count == 3208
    orbit_sizes = Counter(len(members) for members in classes.values())
    assert sum(size * count for size, count in orbit_sizes.items()) == 3208
    rows = []
    class_degrees = Counter()
    for representative, members in sorted(classes.items()):
        orbit = set(d4_images(representative))
        assert orbit == set(members)
        degree = ''.join(map(str, sorted(prior.degrees(representative),
                                         reverse=True)))
        class_degrees[degree] += 1
        degree_two, first_cones = active_profiles(representative)
        rows.append({"representative": [list(edge) for edge in representative],
                     "orbit_size": len(members), "degree_sequence": degree,
                     "degree_two_active_tail_counts": degree_two,
                     "first_commuting_cone_sizes": first_cones})
    result = {"source_survivors": 3208,
              "symmetry": "physical dihedral D4 only",
              "class_count": len(classes),
              "orbit_size_counts": dict(sorted(orbit_sizes.items())),
              "class_degree_counts": dict(sorted(class_degrees.items())),
              "representatives": rows}
    path = ROOT / "global_lower_bound7_dihedral_orbits_result.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("source_survivors", "class_count",
                      "orbit_size_counts", "class_degree_counts")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
