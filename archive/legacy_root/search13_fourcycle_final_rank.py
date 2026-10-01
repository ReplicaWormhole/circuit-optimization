"""Exact cut-rank screen of four-cycle final four-CNOT layers.

Fix the five-CNOT parity prefix family and the exact first 02/13 layer from
runs 198-213. Enumerate the three four-edge cycles on four wires, then all
directed chronological orders (24*16 per cycle). Since cut crossing counts
depend only on the edge multiset, compare exact target Schmidt ranks with
2**crossings before spending effort on local SU(2) fits.
"""

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair, PAIRS, enumerate_near_cases
from search13_crosspair_second_layer_rank import CENTER, four_wire_factor
from search13_crosspair_single_reroute_rank import CUTS, exact_rank, realignment


def cycles():
    all_edges = {tuple(edge) for edge in combinations(range(4), 2)}
    out = []
    for missing in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        cycle = tuple(sorted(all_edges - set(missing)))
        assert len(cycle) == 4
        assert all(sum(wire in edge for edge in cycle) == 2 for wire in range(4))
        out.append({"missing_matching": [list(edge) for edge in missing],
                    "edges": [list(edge) for edge in cycle]})
    return out


def crossing_count(edges, cut):
    return sum((a in cut) != (b in cut) for a, b in edges)


def schedules(cycle):
    edges = [tuple(edge) for edge in cycle["edges"]]
    for order in permutations(edges):
        for flips in product((0, 1), repeat=4):
            yield [edge[::-1] if flip else edge for edge, flip in zip(order, flips)]


def main():
    cases = enumerate_near_cases()
    assert len(cases) == 196
    unique = {}
    for case in cases:
        key = (tuple(case["endpoint"]), case["missing_mask"],
               int(case["missing_coefficient_pi_over_8"]) % 16)
        unique.setdefault(key, 0)
        unique[key] += 1
    assert len(unique) == 51

    exact = CyclotomicPair()
    zero, one, imag = exact.zero, exact.one, exact.imag
    x = [[zero, one], [one, zero]]
    y = [[zero, -imag], [imag, zero]]
    xx, yy = exact.kron(x, x), exact.kron(y, y)
    co, si = exact.trig(Fraction(1, 8))
    f2 = exact.matmul(exact.add_scaled(co, exact.eye4, imag*si, xx),
                      exact.add_scaled(co, exact.eye4, imag*si, yy))
    center = exact.kron(exact.kron(exact.zyz(CENTER[0]), exact.zyz(CENTER[1])),
                        exact.kron(exact.zyz(CENTER[2]), exact.zyz(CENTER[3])))
    baseline_second = exact.matmul(exact.kron(f2, f2), center)

    cycle_records = cycles()
    for cycle in cycle_records:
        cycle["crossings"] = {"".join(map(str, cut)):
                              crossing_count(cycle["edges"], cut) for cut in CUTS}
        cycle["directed_chronological_schedule_count"] = sum(1 for _ in schedules(cycle))
        assert cycle["directed_chronological_schedule_count"] == 384
    rank_histogram = {"".join(map(str, cut)): Counter() for cut in CUTS}
    generator_cache = {}
    targets = []
    compatible = []
    for (endpoint, mask, coefficient), multiplicity in unique.items():
        key = endpoint, mask
        if key not in generator_cache:
            bits = exact.mask_bits(endpoint, mask)
            q02 = exact.transported_factor(endpoint, PAIRS[0], bits)
            q13 = exact.transported_factor(endpoint, PAIRS[1], bits)
            generator_cache[key] = bits, exact.matmul(
                baseline_second, four_wire_factor(q02, q13, exact))
        bits, bq = generator_cache[key]
        ca, sa = exact.trig(Fraction(coefficient, 16))
        target = exact.add_scaled(ca, baseline_second, -imag*sa, bq)
        ranks = {}
        for cut in CUTS:
            label = "".join(map(str, cut))
            ranks[label] = exact_rank(realignment(target, cut, exact), exact)
            rank_histogram[label][ranks[label]] += 1
        matching_cycles = []
        for index, cycle in enumerate(cycle_records):
            if all(ranks[label] <= 2**crossings
                   for label, crossings in cycle["crossings"].items()):
                matching_cycles.append(index)
                compatible.append({"endpoint": list(endpoint),
                                   "missing_mask": mask,
                                   "coefficient_pi_over_8_mod16": coefficient,
                                   "cycle_index": index,
                                   "path_multiplicity": multiplicity})
        targets.append({"endpoint": list(endpoint), "missing_mask": mask,
                        "coefficient_pi_over_8_mod16": coefficient,
                        "post_prefix_wire_mask": bits,
                        "path_multiplicity": multiplicity,
                        "target_schmidt_ranks": ranks,
                        "compatible_cycle_indices": matching_cycles})

    result = {"near_cases": len(cases), "distinct_targets": len(unique),
              "distinct_transported_generators": len(generator_cache),
              "cycle_count": len(cycle_records),
              "directed_chronological_schedule_count": sum(
                  cycle["directed_chronological_schedule_count"] for cycle in cycle_records),
              "cycles": cycle_records,
              "rank_histogram": {key: dict(value) for key, value in rank_histogram.items()},
              "compatible_target_cycle_count": len(compatible),
              "compatible_target_cycles": compatible,
              "targets": targets}
    out = ROOT / "search13_fourcycle_final_rank_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("cycles", "compatible_target_cycles", "targets")},
                     indent=2))


if __name__ == "__main__":
    main()
