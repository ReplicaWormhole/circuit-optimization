"""Exact 03/12 reroute screen for the one-mask-missing five-prefix cases.

The exact first 02/13 layer (including its input frame and pair-local prefix
endpoint map) is fixed. After transporting the omitted parity phase through
that layer, ask whether the *entire* target second-layer operator can factor
as one two-qubit unitary on 03 times one on 12. Such a factor has operator
Schmidt rank one across 03|12. Exact nonzero realignment minors exclude it.
"""

import json
from collections import Counter
from fractions import Fraction

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair, PAIRS, enumerate_near_cases


CENTER = {
    0: (Fraction(1, 2), Fraction(0), Fraction(-1)),
    1: (Fraction(1), Fraction(0), Fraction(3, 2)),
    2: (Fraction(1), Fraction(0), Fraction(3, 2)),
    3: (Fraction(1, 2), Fraction(-3, 2), Fraction(0)),
}


def four_wire_factor(q02, q13, exact):
    out = [[exact.zero for _ in range(16)] for _ in range(16)]
    for y in range(16):
        y0, y1, y2, y3 = ((y >> shift) & 1 for shift in (3, 2, 1, 0))
        out02, out13 = (y0 << 1) | y2, (y1 << 1) | y3
        for x in range(16):
            x0, x1, x2, x3 = ((x >> shift) & 1 for shift in (3, 2, 1, 0))
            in02, in13 = (x0 << 1) | x2, (x1 << 1) | x3
            out[y][x] = q02[out02][in02] * q13[out13][in13]
    return out


def rank_one_across_03_12(operator, exact):
    realigned = [[exact.zero for _ in range(16)] for _ in range(16)]
    for y in range(16):
        out03 = (y & 8) >> 2 | (y & 1)
        out12 = (y & 4) >> 1 | (y & 2) >> 1
        for x in range(16):
            in03 = (x & 8) >> 2 | (x & 1)
            in12 = (x & 4) >> 1 | (x & 2) >> 1
            realigned[4*out03+in03][4*out12+in12] = operator[y][x]
    pivot = next(((r, c) for r in range(16) for c in range(16)
                  if realigned[r][c] != exact.zero), None)
    assert pivot is not None
    r0, c0 = pivot
    minor = next(((r, c) for r in range(16) for c in range(16)
                  if realigned[r0][c0]*realigned[r][c]
                  - realigned[r0][c]*realigned[r][c0] != exact.zero), None)
    return minor is None, {"pivot": [r0, c0],
                           "nonzero_minor_rows_columns":
                               [[r0, minor[0]], [c0, minor[1]]] if minor else None}


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
    ca, sa = exact.trig(Fraction(1, 8))
    f2 = exact.matmul(exact.add_scaled(ca, exact.eye4, imag*sa, xx),
                      exact.add_scaled(ca, exact.eye4, imag*sa, yy))
    last_layer = exact.kron(f2, f2)
    center = exact.kron(exact.kron(exact.zyz(CENTER[0]), exact.zyz(CENTER[1])),
                        exact.kron(exact.zyz(CENTER[2]), exact.zyz(CENTER[3])))
    baseline_second = exact.matmul(last_layer, center)
    eye16 = exact.identity(16)

    phase_cache = {}
    records = []
    phase_rank_one, full_rank_one = [], []
    for (endpoint, mask, coefficient), multiplicity in unique.items():
        key = endpoint, mask
        if key not in phase_cache:
            bits = exact.mask_bits(endpoint, mask)
            q02 = exact.transported_factor(endpoint, PAIRS[0], bits)
            q13 = exact.transported_factor(endpoint, PAIRS[1], bits)
            q16 = four_wire_factor(q02, q13, exact)
            bq = exact.matmul(baseline_second, q16)
            phase_cache[key] = bits, q16, bq
        bits, q16, bq = phase_cache[key]
        co, si = exact.trig(Fraction(coefficient, 16))
        phase = exact.add_scaled(co, eye16, -imag*si, q16)
        full = exact.add_scaled(co, baseline_second, -imag*si, bq)
        phase_factors, phase_witness = rank_one_across_03_12(phase, exact)
        full_factors, full_witness = rank_one_across_03_12(full, exact)
        record = {"endpoint": list(endpoint), "missing_mask": mask,
                  "coefficient_pi_over_8_mod16": coefficient,
                  "post_prefix_wire_mask": bits, "multiplicity": multiplicity,
                  "transported_phase_rank_one_03_12": phase_factors,
                  "full_second_layer_rank_one_03_12": full_factors,
                  "phase_rank_witness": phase_witness,
                  "full_rank_witness": full_witness}
        records.append(record)
        if phase_factors:
            phase_rank_one.append(record)
        if full_factors:
            full_rank_one.append(record)

    result = {"near_cases": len(cases), "distinct_phase_operators": len(unique),
              "distinct_transported_generators": len(phase_cache),
              "phase_rank_one_count": len(phase_rank_one),
              "full_second_layer_rank_one_count": len(full_rank_one),
              "phase_rank_one_cases": phase_rank_one,
              "full_second_layer_rank_one_cases": full_rank_one,
              "records": records}
    out = ROOT / "search13_crosspair_second_layer_rank_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "records"}, indent=2))


if __name__ == "__main__":
    main()
