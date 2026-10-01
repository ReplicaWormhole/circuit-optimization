"""Exact mixed-cut screen of the three untested overlapping block pairs.

The two first-layer blocks and two final-layer blocks each commute. For each
first/final overlap, the relevant three-wire factor is F_final L_shared F_first;
all other center locals move to its exterior and leave its CNOT cost unchanged.
This is a fixed-operator necessary test, not a global circuit lower bound.
"""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path

from search13_middle_threequbit import (
    PAIRS, MASKS, PRIMES, embed_pair, graph_edges, mincut,
    reduce_field, realign, rank_mod,
)
from search13_secondpair_mask_transport import CyclotomicPair
import sympy as sp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_sibling_middle_blocks_result.json"
CENTER = {
    0: (Fraction(1, 2), Fraction(0), Fraction(-1)),
    1: (Fraction(1), Fraction(0), Fraction(3, 2)),
    2: (Fraction(1), Fraction(0), Fraction(3, 2)),
    3: (Fraction(1, 2), Fraction(-3, 2), Fraction(0)),
}
CASES = (("F01_L0_F02", (0, 1), (0, 2), (0, 1, 2)),
         ("F23_L2_F02", (2, 3), (0, 2), (0, 2, 3)),
         ("F23_L3_F13", (2, 3), (1, 3), (1, 2, 3)))


def final_pair(exact):
    c, s = exact.trig(Fraction(1, 8))
    x = [[exact.zero, exact.one], [exact.one, exact.zero]]
    y = [[exact.zero, -exact.imag], [exact.imag, exact.zero]]
    return exact.matmul(exact.add_scaled(c, exact.eye4, exact.imag*s,
                                         exact.kron(x, x)),
                        exact.add_scaled(c, exact.eye4, exact.imag*s,
                                         exact.kron(y, y)))


def target_for_case(exact, final, first, wires):
    local = exact.identity(1)
    for wire in wires:
        factor = exact.zyz(CENTER[wire]) if wire in set(final) & set(first) else exact.eye2
        local = exact.kron(local, factor)
    positions = lambda pair: tuple(wires.index(w) for w in pair)
    return exact.matmul(embed_pair(exact, final_pair(exact), positions(final)),
                        exact.matmul(local, embed_pair(exact, exact.f,
                                                       positions(first))))


def modular_ranks(target):
    ranks = {mask: 0 for mask in MASKS}
    by_prime = {}
    for prime in PRIMES:
        root = pow(int(sp.primitive_root(prime)), (prime - 1)//96, prime)
        assert (pow(root, 32, prime) - pow(root, 16, prime) + 1) % prime == 0
        matrix = [[reduce_field(v, root, prime) for v in row] for row in target]
        by_prime[str(prime)] = {str(mask): rank_mod(realign(matrix, mask), prime)
                                for mask in MASKS}
        for mask in MASKS:
            ranks[mask] = max(ranks[mask], by_prime[str(prime)][str(mask)])
    return ranks, by_prime


def main():
    exact = CyclotomicPair()
    results = []
    for name, final, first, wires in CASES:
        target = target_for_case(exact, final, first, wires)
        ranks, by_prime = modular_ranks(target)
        exclusions = []
        survivors = []
        for schedule in product(PAIRS, repeat=3):
            edges = graph_edges(schedule)
            witness = next(((mask, ranks[mask], mincut(edges, mask))
                            for mask in MASKS
                            if ranks[mask] > 2**mincut(edges, mask)), None)
            if witness is None:
                survivors.append([list(pair) for pair in schedule])
            else:
                exclusions.append({"schedule": [list(pair) for pair in schedule],
                                   "witness": witness})
        results.append({"name": name, "wires": wires,
                        "rank_lower_bounds": {str(k): v for k, v in ranks.items()},
                        "ranks_by_prime": by_prime,
                        "excluded_count": len(exclusions), "excluded": exclusions,
                        "surviving_schedules": survivors})
        print(name, "excluded", len(exclusions), "survive", len(survivors), flush=True)
    output = {"field": "Q(zeta_96)", "primes": PRIMES,
              "boundary_order": "out0,out1,out2,in0,in1,in2 on each case's wires",
              "scope": "fixed overlapping three-wire blocks only; changed intermediate gauges and global 13-CNOT architectures remain open",
              "cases": results}
    OUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
