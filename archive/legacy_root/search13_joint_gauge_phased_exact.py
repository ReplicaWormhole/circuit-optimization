"""Exact cut-rank certificates for selected complex double-row gauges.

The matrices and ranks are computed in Q(zeta_96).  These certify ranks of
two fixed gauged tail operators, not a lower bound for arbitrary gauges.
"""

import json
from fractions import Fraction

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair
from search13_crosspair_second_layer_rank import CENTER, four_wire_factor
from search13_crosspair_single_reroute_rank import CUTS, exact_rank, realignment


PAIRS = ((4, 7), (8, 11))
SPECS = (("best_full", (Fraction(-1, 2), Fraction(-1, 2))),
         ("best_final", (Fraction(-1, 4), Fraction(1, 4))))


def base_matrices(exact):
    zero, one, imag = exact.zero, exact.one, exact.imag
    x = [[zero, one], [one, zero]]
    y = [[zero, -imag], [imag, zero]]
    xx, yy = exact.kron(x, x), exact.kron(y, y)
    co, si = exact.trig(Fraction(1, 8))
    f2 = exact.matmul(exact.add_scaled(co, exact.eye4, imag*si, xx),
                      exact.add_scaled(co, exact.eye4, imag*si, yy))
    center = exact.kron(exact.kron(exact.zyz(CENTER[0]), exact.zyz(CENTER[1])),
                        exact.kron(exact.zyz(CENTER[2]), exact.zyz(CENTER[3])))
    final = exact.matmul(exact.kron(f2, f2), center)
    first = four_wire_factor(exact.f, exact.f, exact)
    input_local = exact.kron(
        exact.kron(exact.local[0], exact.local[1]),
        exact.kron(exact.local[2], exact.local[3]))
    full = exact.matmul(exact.matmul(final, first), input_local)
    return full, final


def gauged(matrix, phases, exact):
    out = [row[:] for row in matrix]
    co, si = exact.trig(Fraction(1, 4))
    for (first, second), phase in zip(PAIRS, phases):
        exponent = 48*phase
        assert exponent.denominator == 1
        ep = exact.z**int(exponent)
        a, b = out[first][:], out[second][:]
        out[first] = [co*u-ep**-1*si*v for u, v in zip(a, b)]
        out[second] = [ep*si*u+co*v for u, v in zip(a, b)]
    return out


def ranks(matrix, exact):
    return {"".join(map(str, cut)):
            exact_rank(realignment(matrix, cut, exact), exact) for cut in CUTS}


def main():
    exact = CyclotomicPair()
    full, final = base_matrices(exact)
    result = {"field": "Q(zeta_96)", "pairs": PAIRS,
              "rotation_angle_pi_multiple": "1/4", "rows": []}
    for label, phases in SPECS:
        gf, gl = gauged(full, phases, exact), gauged(final, phases, exact)
        row = {"name": label,
               "phase_pi_multiples": [str(p) for p in phases],
               "full_ranks": ranks(gf, exact),
               "final_ranks": ranks(gl, exact)}
        result["rows"].append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    path = ROOT / "search13_joint_gauge_phased_exact_result.json"
    path.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
