"""Exact operator-Schmidt-rank obstruction for second-pair phase uptake.

For each distinct one-mask-missing case, form the omitted phase after the
fixed first 02/13 layer in Q(zeta_96). An independent 01 and 23 two-qubit
layer has Schmidt rank one across 01|23. Exhibit a nonzero 2x2 realignment
minor for every transported phase to exclude such a factorization.
"""

import json
from collections import Counter
from fractions import Fraction

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair, PAIRS, enumerate_near_cases


def main():
    cases = enumerate_near_cases()
    assert len(cases) == 196
    unique = {}
    for case in cases:
        key = (tuple(case["endpoint"]), case["missing_mask"],
               int(case["missing_coefficient_pi_over_8"]) % 16)
        unique.setdefault(key, 0)
        unique[key] += 1

    exact = CyclotomicPair()
    zero, one, imag = exact.zero, exact.one, exact.imag
    results = []
    rank_one = []
    for (endpoint, mask, coefficient), multiplicity in unique.items():
        bits = exact.mask_bits(endpoint, mask)
        q02 = exact.transported_factor(endpoint, PAIRS[0], bits)
        q13 = exact.transported_factor(endpoint, PAIRS[1], bits)
        co, si = exact.trig(Fraction(coefficient, 16))
        realigned = [[zero for _ in range(16)] for _ in range(16)]
        for output in range(16):
            y0, y1, y2, y3 = ((output >> s) & 1 for s in (3, 2, 1, 0))
            out02, out13 = (y0 << 1) | y2, (y1 << 1) | y3
            for input_index in range(16):
                x0, x1, x2, x3 = ((input_index >> s) & 1 for s in (3, 2, 1, 0))
                in02, in13 = (x0 << 1) | x2, (x1 << 1) | x3
                value = (co * (one if output == input_index else zero)
                         - imag * si * q02[out02][in02] * q13[out13][in13])
                row = (output >> 2)*4 + (input_index >> 2)
                col = (output & 3)*4 + (input_index & 3)
                realigned[row][col] = value

        pivot = next(((r, c) for r in range(16) for c in range(16)
                      if realigned[r][c] != zero), None)
        assert pivot is not None
        r0, c0 = pivot
        witness = next(((r, c) for r in range(16) for c in range(16)
                        if realigned[r0][c0]*realigned[r][c]
                        - realigned[r0][c]*realigned[r][c0] != zero), None)
        row = {"endpoint": list(endpoint), "missing_mask": mask,
               "coefficient_pi_over_8_mod16": coefficient,
               "multiplicity": multiplicity,
               "post_prefix_wire_mask": bits,
               "nonzero_minor_rows_columns":
                   [[r0, witness[0]], [c0, witness[1]]] if witness else None}
        results.append(row)
        if witness is None:
            rank_one.append(row)

    result = {"case_count": len(cases), "distinct_phase_operators": len(unique),
              "coefficient_histogram": dict(sorted(Counter(
                  int(case["missing_coefficient_pi_over_8"]) % 16
                  for case in cases).items())),
              "schmidt_rank_at_least_two": len(results)-len(rank_one),
              "rank_one_cases": rank_one, "results": results}
    out = ROOT / "search13_secondpair_phase_rank_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "results"}, indent=2))


if __name__ == "__main__":
    main()
