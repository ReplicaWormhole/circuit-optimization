"""Exact second-pair transport screen for one omitted parity phase.

Enumerate five-CNOT parity tours with exactly one omitted required mask and
a first-pair-local linear mismatch. Conjugate that mask through the fixed
first matchgate layer F02 F13 and its exact input local frame, including the
linear mismatch. Test exactly over Q(zeta_96) whether it becomes supported
on 01 or 23 before the center local layer. Local center gates preserve wire
support, so failure is an exact necessary obstruction to second-pair uptake.
"""

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import sympy as sp

from search13_firstpair_mask_relocation import (
    BASIS, EXACT6_ROWS, ROOT, pair_local_endpoints,
)


PAIRS = ((0, 2), (1, 3))
SECOND_PAIRS = ((0, 1), (2, 3))
LOCAL = {
    0: (Fraction(3, 4), Fraction(0), Fraction(0)),
    1: (Fraction(1, 2), Fraction(-1, 2), Fraction(-1, 3)),
    2: (Fraction(1, 2), Fraction(-3, 2), Fraction(1, 6)),
    3: (Fraction(3, 4), Fraction(-1), Fraction(1)),
}


def enumerate_near_cases():
    models = []
    for family in ("fourmask", "fivemask"):
        data = json.loads((ROOT / f"algebraic14_{family}_models.json").read_text())
        for index, model in enumerate(data):
            models.append((family, index, frozenset(model["support"]),
                           model["coefficients_pi_over_8"]))
    targets = pair_local_endpoints()
    cases = []

    def walk(rows, seen, path):
        if len(path) == 5:
            for active_pair in targets.get(rows, ()):
                for family, index, support, coeffs in models:
                    missing = support - seen
                    if len(missing) != 1:
                        continue
                    mask = next(iter(missing))
                    cases.append({"family": family, "model_index": index,
                                  "support": sorted(support),
                                  "missing_mask": mask,
                                  "missing_coefficient_pi_over_8": coeffs[str(mask)],
                                  "active_pair": list(active_pair),
                                  "endpoint": list(rows),
                                  "path": [list(edge) for edge in path]})
            return
        for a in range(4):
            for b in range(4):
                if a == b:
                    continue
                next_rows = list(rows)
                next_rows[b] ^= rows[a]
                mask = next_rows[b]
                next_seen = seen | {mask} if mask.bit_count() > 1 else seen
                walk(tuple(next_rows), next_seen, path + ((a, b),))

    walk(BASIS, frozenset(), ())
    return cases


class CyclotomicPair:
    def __init__(self):
        self.field = sp.QQ.cyclotomic_field(96)
        self.z = self.field.unit
        self.zero = self.field.zero
        self.one = self.field.one
        self.half = self.field.from_sympy(sp.Rational(1, 2))
        self.imag = self.z**24
        self.eye2 = self.identity(2)
        self.eye4 = self.identity(4)
        self.zz = [[self.one, self.zero], [self.zero, -self.one]]
        x = [[self.zero, self.one], [self.one, self.zero]]
        y = [[self.zero, -self.imag], [self.imag, self.zero]]
        xx, yy = self.kron(x, x), self.kron(y, y)
        ca, sa = self.trig(Fraction(1, 4))
        cb, sb = self.trig(Fraction(1, 8))
        self.f = self.matmul(self.add_scaled(ca, self.eye4, self.imag*sa, xx),
                             self.add_scaled(cb, self.eye4, self.imag*sb, yy))
        self.fdag = self.matmul(self.add_scaled(cb, self.eye4, -self.imag*sb, yy),
                                self.add_scaled(ca, self.eye4, -self.imag*sa, xx))
        self.local = {wire: self.zyz(angles) for wire, angles in LOCAL.items()}
        self.local_dag = {wire: self.zyz_dag(angles) for wire, angles in LOCAL.items()}

    def identity(self, n):
        return [[self.one if i == j else self.zero for j in range(n)]
                for i in range(n)]

    def matmul(self, a, b):
        return [[sum((a[i][k]*b[k][j] for k in range(len(b))), self.zero)
                 for j in range(len(b[0]))] for i in range(len(a))]

    @staticmethod
    def transpose(a):
        return [list(row) for row in zip(*a)]

    def kron(self, a, b):
        na, nb = len(a), len(b)
        return [[a[i//nb][j//nb]*b[i%nb][j%nb]
                 for j in range(na*nb)] for i in range(na*nb)]

    def add_scaled(self, a, x, b, y):
        return [[a*x[i][j]+b*y[i][j] for j in range(len(x))]
                for i in range(len(x))]

    def trig(self, theta):
        k = 48*theta
        assert k.denominator == 1
        p, n = self.z**int(k), self.z**-int(k)
        return (p+n)*self.half, (p-n)*self.half/self.imag

    def rz(self, theta):
        k = 24*theta
        assert k.denominator == 1
        return [[self.z**-int(k), self.zero], [self.zero, self.z**int(k)]]

    def ry(self, theta):
        k = 24*theta
        assert k.denominator == 1
        p, n = self.z**int(k), self.z**-int(k)
        c, s = (p+n)*self.half, (p-n)*self.half/self.imag
        return [[c, -s], [s, c]]

    def zyz(self, angles):
        theta, phi, lam = angles
        return self.matmul(self.matmul(self.rz(phi), self.ry(theta)), self.rz(lam))

    def zyz_dag(self, angles):
        theta, phi, lam = angles
        return self.matmul(self.matmul(self.rz(-lam), self.ry(-theta)), self.rz(-phi))

    def permutation(self, endpoint, pair):
        """Two-wire P with E6_pair=P E5_pair as binary row maps."""
        a, b = pair
        e0, e1 = endpoint[a], endpoint[b]

        def coefficients(mask):
            for c0 in (0, 1):
                for c1 in (0, 1):
                    if (e0 if c0 else 0) ^ (e1 if c1 else 0) == mask:
                        return c0, c1
            raise AssertionError("endpoint is not pair-local")

        upper = coefficients(EXACT6_ROWS[a])
        lower = coefficients(EXACT6_ROWS[b])
        p = [[self.zero for _ in range(4)] for _ in range(4)]
        for col in range(4):
            x0, x1 = col >> 1, col & 1
            y0 = (upper[0]*x0) ^ (upper[1]*x1)
            y1 = (lower[0]*x0) ^ (lower[1]*x1)
            p[(y0 << 1) | y1][col] = self.one
        return p

    def mask_bits(self, endpoint, mask):
        for bits in range(16):
            total = 0
            for wire in range(4):
                if bits & (8 >> wire):
                    total ^= endpoint[wire]
            if total == mask:
                return bits
        raise AssertionError("linear endpoint is not invertible")

    def transported_factor(self, endpoint, pair, bits):
        a, b = pair
        pauli = self.kron(self.zz if bits & (8 >> a) else self.eye2,
                          self.zz if bits & (8 >> b) else self.eye2)
        p = self.permutation(endpoint, pair)
        l = self.kron(self.local[a], self.local[b])
        ldag = self.kron(self.local_dag[a], self.local_dag[b])
        u = self.matmul(self.matmul(self.f, l), p)
        udag = self.matmul(self.matmul(self.transpose(p), ldag), self.fdag)
        return self.matmul(self.matmul(u, pauli), udag)

    def local_support(self, q):
        """Exact possible single-wire support positions within a pair."""
        o0 = [[(q[2*i][2*j]+q[2*i+1][2*j+1])*self.half
               for j in range(2)] for i in range(2)]
        o1 = [[(q[i][j]+q[i+2][j+2])*self.half
               for j in range(2)] for i in range(2)]
        possible = []
        if q == self.kron(o0, self.eye2):
            possible.append(0)
        if q == self.kron(self.eye2, o1):
            possible.append(1)
        return possible


def main():
    cases = enumerate_near_cases()
    assert len(cases) == 196
    exact = CyclotomicPair()
    cache = {}
    survivors = []
    support_histogram = Counter()
    for case in cases:
        endpoint = tuple(case["endpoint"])
        bits = exact.mask_bits(endpoint, case["missing_mask"])
        profile = {}
        for pair in PAIRS:
            key = (endpoint, pair, bits & sum(8 >> wire for wire in pair))
            if key not in cache:
                q = exact.transported_factor(endpoint, pair, bits)
                cache[key] = exact.local_support(q)
            profile[pair] = cache[key]
        support_histogram[str((tuple(profile[PAIRS[0]]), tuple(profile[PAIRS[1]])))] += 1
        target_pairs = []
        if 0 in profile[PAIRS[0]] and 0 in profile[PAIRS[1]]:
            target_pairs.append((0, 1))
        if 1 in profile[PAIRS[0]] and 1 in profile[PAIRS[1]]:
            target_pairs.append((2, 3))
        if target_pairs:
            survivors.append({**case, "post_prefix_wire_mask": bits,
                              "second_pairs": [list(pair) for pair in target_pairs]})
    result = {"ordered_five_cnot_paths": 12**5, "near_case_count": len(cases),
              "distinct_exact_factor_checks": len(cache),
              "support_profile_histogram": dict(support_histogram),
              "second_pair_survivor_count": len(survivors),
              "survivors": survivors}
    out = ROOT / "search13_secondpair_mask_transport_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "survivors"}, indent=2))


if __name__ == "__main__":
    main()
