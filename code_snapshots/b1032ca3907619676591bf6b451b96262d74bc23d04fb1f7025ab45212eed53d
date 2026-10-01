"""Exact modular commutant screen for a seven-CNOT tail after the fixed prefix.

For a tail T diagonalizing A=P V4 P^{-1}, every T^{-1} Z_j T commutes
with A. A forward support cone in the reversed tail must therefore contain
a nonzero traceless commuting operator. All ranks are computed over two
prime fields containing a primitive sixteenth root; full column rank at
either prime certifies full rank in characteristic zero.
"""

from collections import Counter
from itertools import combinations, product
import json
from math import pi
from pathlib import Path

from sympy import primitive_root

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_fulltail_invariant_commutant_result.json"
EDGES = tuple(combinations(range(4), 2))
PRIMES = (97, 193)


def prefix_action():
    gates = json.loads(SOURCE.read_text())["gates"]
    perm = list(range(16))
    exponent = [0] * 16
    count = 0
    for gate in gates:
        if gate["gate"] == "cx":
            count += 1
            mask = 1 << (3 - gate["target"])
            test = 1 << (3 - gate["control"])
            perm = [y ^ mask if y & test else y for y in perm]
        elif gate["gate"] == "rz":
            # Remove an irrelevant scalar: Rz(k*pi/8) ~ diag(1,zeta16^k).
            from check_circuit import angle
            raw = angle(gate["theta"]) * 8 / pi
            k = round(raw)
            assert abs(raw - k) < 1e-12
            test = 1 << (3 - gate["qubit"])
            exponent = [(e + k * bool(y & test)) % 16
                        for e, y in zip(exponent, perm)]
        else:
            raise ValueError("unexpected nonmonomial prefix gate")
        if count == 6:
            break
    assert count == 6 and len(set(perm)) == 16
    return perm, exponent


def residual_shift_action(prefix_perm, prefix_exp):
    inverse = [0] * 16
    for x, y in enumerate(prefix_perm):
        inverse[y] = x
    shift = lambda x: ((x & 1) << 3) | (x >> 1)
    perm = [prefix_perm[shift(inverse[y])] for y in range(16)]
    phase = [(prefix_exp[shift(inverse[y])] - prefix_exp[inverse[y]]) % 16
             for y in range(16)]
    return perm, phase


def pauli_action(label, imaginary, prime):
    flip = sum(1 << (3 - j) for j, a in enumerate(label) if a in "XY")
    coefficients = []
    for x in range(16):
        coeff = 1
        for j, a in enumerate(label):
            bit = (x >> (3 - j)) & 1
            if a == "Z" and bit:
                coeff = -coeff
            if a == "Y":
                coeff *= -imaginary if bit else imaginary
        coefficients.append(coeff % prime)
    return flip, coefficients


def rank_mod(columns, prime):
    if not columns:
        return 0
    rows = [list(row) for row in zip(*columns)]
    rank = 0
    for col in range(len(columns)):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        factor = pow(rows[rank][col], -1, prime)
        rows[rank] = [(v * factor) % prime for v in rows[rank]]
        for i in range(rank + 1, len(rows)):
            if rows[i][col]:
                factor = rows[i][col]
                rows[i] = [(a - factor * b) % prime
                           for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def commutant_ranks(res_perm, phase_exp, subset, prime):
    root = pow(primitive_root(prime), (prime - 1) // 16, prime)
    assert pow(root, 16, prime) == 1 and pow(root, 8, prime) != 1
    imag = pow(root, 4, prime)
    phase = [pow(root, e, prime) for e in phase_exp]
    labels = ["".join(axes) for axes in product("IXYZ", repeat=4)
              if any(a != "I" for a in axes)
              and all(a == "I" for j, a in enumerate(axes) if j not in subset)]
    columns = []
    for label in labels:
        flip, coeff = pauli_action(label, imag, prime)
        comm = [0] * 256
        for x in range(16):
            y = x ^ flip
            comm[16 * res_perm[y] + x] += phase[y] * coeff[x]
            z = res_perm[x]
            comm[16 * (z ^ flip) + x] -= phase[x] * coeff[z]
        columns.append([v % prime for v in comm])
    return len(labels), rank_mod(columns, prime)


def cone(schedule, wire):
    current = {wire}
    for edge in reversed(schedule):
        if current.intersection(edge):
            current.update(edge)
    return tuple(sorted(current))


def main():
    prefix = prefix_action()
    residual = residual_shift_action(*prefix)
    ranks = {}
    for size in range(1, 4):
        for subset in combinations(range(4), size):
            entries = [commutant_ranks(*residual, subset, prime)
                       for prime in PRIMES]
            ranks["".join(map(str, subset))] = {
                "basis_dimension": entries[0][0],
                "ranks_mod_prime": {str(p): r for p, (_, r) in zip(PRIMES, entries)},
            }
    # A full column rank at either prime rigorously excludes a commuting
    # traceless Hermitian supported in that subset.
    forbidden = {tuple(map(int, key)) for key, row in ranks.items()
                 if max(row["ranks_mod_prime"].values()) == row["basis_dimension"]}
    counts = Counter()
    excluded_examples = {}
    for schedule in product(EDGES, repeat=7):
        witnesses = [j for j in range(4) if cone(schedule, j) in forbidden]
        if witnesses:
            counts["excluded_by_commutant_cone"] += 1
            excluded_examples.setdefault(str(len(witnesses)), {
                "schedule": schedule, "witnesses": witnesses})
        else:
            counts["survives_commutant_cone"] += 1
    assert sum(counts.values()) == 6**7
    result = {
        "fixed_prefix_source": SOURCE.name,
        "prefix_cnot_count": 6,
        "field_primes": PRIMES,
        "subset_commutant_ranks": ranks,
        "certified_no_commutant_subsets": [list(s) for s in sorted(forbidden)],
        "seven_cnot_unoriented_schedules": 6**7,
        "exclusive_counts": dict(counts),
        "excluded_examples": excluded_examples,
        "scope": "necessary support-cone test only; no seven-CNOT synthesis or no-go",
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"counts": counts, "forbidden": sorted(forbidden)}, default=list))


if __name__ == "__main__":
    main()
