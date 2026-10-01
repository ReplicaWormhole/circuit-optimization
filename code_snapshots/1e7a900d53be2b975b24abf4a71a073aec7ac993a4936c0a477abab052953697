"""Exact integer-valued spectral-path scan for V4 Pauli products.

All matrices have Gaussian-integer entries; products are small enough that
complex128 addition/multiplication is exact. Projectors are represented as 4P.
This is a finite diagnostic, not a classification of arbitrary local axes.
"""

import json
from itertools import permutations, product

import numpy as np


PAULI = {
    "I": np.array([[1, 0], [0, 1]], dtype=np.complex128),
    "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
}


def kron_all(items):
    out = np.array([[1]], dtype=np.complex128)
    for item in items:
        out = np.kron(out, item)
    return out


def main():
    shift = np.zeros((16, 16), dtype=np.complex128)
    for x in range(16):
        shift[((x & 1) << 3) | (x >> 1), x] = 1
    powers = [np.linalg.matrix_power(shift, r) for r in range(4)]
    projector4 = [sum(((-1j) ** (k * r) * powers[r] for r in range(4)),
                      np.zeros((16, 16), dtype=np.complex128)) for k in range(4)]
    assert [int(np.trace(p).real / 4) for p in projector4] == [6, 3, 4, 3]
    triples = list(permutations(range(4), 3))
    result = {"normalization": "Mk=4Pk", "single": {}, "pairs": {}}
    for support_size in (1, 2):
        from itertools import combinations
        for support in combinations(range(4), support_size):
            passing = []
            failing = {}
            for letters in product("XYZ", repeat=support_size):
                label = "".join(f"{a}{j}" for j, a in zip(support, letters))
                factors = [PAULI[letters[support.index(j)]] if j in support else PAULI["I"]
                           for j in range(4)]
                op = kron_all(factors)
                defects = [projector4[a] @ op @ projector4[b] @ op @ projector4[c]
                           for a, b, c in triples]
                max_abs = int(max(np.max(np.abs(d)) for d in defects))
                if max_abs == 0:
                    passing.append(label)
                else:
                    failing[label] = max_abs
            key = "".join(map(str, support))
            result["single" if support_size == 1 else "pairs"][key] = {
                "passing": passing,
                "failing_count": len(failing),
                "minimum_nonzero_defect": min(failing.values()) if failing else None,
            }
    with open("spectral_path4_pair_pauli_scan_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
