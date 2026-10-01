"""Exhaust local Clifford frames after deleting exact15's middle CZ.

Hold all four exact two-qubit pair blocks fixed. Between their two disjoint
pair layers, scan one projective single-qubit Clifford on each of the four
wires (24^4 choices). This tests a noncommuting local-frame replacement of
the CZ, distinct from the mid-layer Z-phase grid. It is a restricted family.
"""

import itertools
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate, right_shift
from input_uniform_clifford14 import H, S, clifford_words
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"
BEST = ROOT / "exact14_mid_clifford_grid_best.json"
RESULT = ROOT / "exact14_mid_clifford_grid_result.json"


def word_matrix(word):
    out = np.eye(2, dtype=complex)
    for letter in word:
        out = (H if letter == "h" else S) @ out
    return out


def main():
    source = json.loads(SOURCE.read_text())
    gates = source["gates"]
    assert gates[49:52] == [{"gate": "h", "qubit": 1},
                            {"gate": "cx", "control": 2, "target": 1},
                            {"gate": "h", "qubit": 1}]
    before = circuit_matrix(gates[:49])
    after = circuit_matrix(gates[52:])
    W = before @ right_shift(4) @ before.conj().T
    words = clifford_words()
    matrices = [word_matrix(word) for word in words]
    pairs = [(i, j, np.kron(matrices[i], matrices[j]))
             for i, j in itertools.product(range(24), repeat=2)]
    best_error = float("inf")
    best_indices = None
    trials = 0
    for i, j, first in pairs:
        for k, l, second in pairs:
            mid = np.kron(first, second)
            conjugated = mid @ W @ mid.conj().T
            rotated = after @ conjugated @ after.conj().T
            error = float(np.max(np.abs(rotated - np.diag(np.diag(rotated)))))
            trials += 1
            if error < best_error:
                best_error = error
                best_indices = (i, j, k, l)
    inserted = []
    for q, index in enumerate(best_indices):
        inserted.extend({"gate": name, "qubit": q} for name in words[index])
    candidate = {"n": 4, "gates": gates[:49] + inserted + gates[52:]}
    BEST.write_text(json.dumps(candidate, indent=2) + "\n")
    checked = evaluate(candidate)
    result = {"trials": trials, "best_indices": best_indices,
              "best_words": ["".join(words[i]) for i in best_indices],
              "best_error": best_error,
              "checked_error": checked["off_diagonal_error"],
              "valid": checked["valid_diagonalizer"],
              "cnot_count": checked["cnot_count"]}
    RESULT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
