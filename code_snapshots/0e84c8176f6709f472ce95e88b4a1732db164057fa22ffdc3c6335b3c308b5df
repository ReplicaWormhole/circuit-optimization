"""Bounded Qiskit scan of global one-qubit Clifford input symmetries.

For any one-qubit C, C^tensor4 commutes with the cyclic wire shift. Thus
U C^tensor4 is another exact diagonalizer whenever U is. This scans the 24
projective one-qubit Cliffords using one Qiskit pass per representative.
The scan is heuristic compilation, not an optimality claim.
"""

import json
from collections import deque
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"
BEST = ROOT / "input_uniform_clifford14_best.json"
RESULT = ROOT / "input_uniform_clifford14_result.json"
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)


def projective_key(matrix):
    vector = matrix.ravel()
    pivot = next(value for value in vector if abs(value) > 1e-9)
    normalized = matrix * np.conjugate(pivot) / abs(pivot)
    return tuple(np.round(normalized.real.ravel(), 9)) + tuple(
        np.round(normalized.imag.ravel(), 9))


def clifford_words():
    seen = {projective_key(np.eye(2, dtype=complex))}
    queue = deque([((), np.eye(2, dtype=complex))])
    words = []
    while queue:
        word, matrix = queue.popleft()
        words.append(word)
        for letter, gate in (("h", H), ("s", S)):
            product = gate @ matrix
            key = projective_key(product)
            if key not in seen:
                seen.add(key)
                queue.append((word + (letter,), product))
    assert len(words) == 24, len(words)
    return words


def main():
    incumbent = to_qiskit(json.loads(SOURCE.read_text()))
    rows = []
    best = None
    for word in clifford_words():
        input_circuit = QuantumCircuit(4)
        for letter in word:
            for qubit in range(4):
                getattr(input_circuit, letter)(qubit)
        source = input_circuit.compose(incumbent)
        optimized = transpile(source, basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        candidate = from_qiskit(optimized)
        check = evaluate(candidate)
        row = {"word": "".join(word), "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid": check["valid_diagonalizer"]}
        rows.append(row)
        if check["valid_diagonalizer"] and (
                best is None or row["cnot_count"] < best["cnot_count"]):
            best = row
            BEST.write_text(json.dumps(candidate, indent=2) + "\n")
    RESULT.write_text(json.dumps({"trials": len(rows), "best": best,
                                  "rows": rows}, indent=2) + "\n")
    print(json.dumps({"trials": len(rows), "best": best}, indent=2))


if __name__ == "__main__":
    main()
