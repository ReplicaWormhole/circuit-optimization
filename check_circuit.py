"""Check a proposed ancilla-free CNOT/one-qubit cycle diagonalizer.

Gates in the input JSON are listed in chronological order. Qubit 0 is the
most-significant bit, matching Equation (1.1) of the cycle paper.
"""

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np


ANGLE = re.compile(r"^([+-]?)(?:(\d+(?:\.\d+)?)\*?)?pi(?:/(\d+))?$")


def angle(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    if not isinstance(value, str):
        raise ValueError(f"invalid angle: {value!r}")
    match = ANGLE.fullmatch(value.replace(" ", ""))
    if match is None:
        raise ValueError(f"invalid angle: {value!r}; use a number or e.g. '3*pi/4'")
    sign, factor, denominator = match.groups()
    return (-1 if sign == "-" else 1) * float(factor or 1) * math.pi / int(denominator or 1)


def one_qubit_matrix(gate):
    name = gate["gate"].lower()
    if name == "h":
        return np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
    if name == "x":
        return np.array([[0, 1], [1, 0]], dtype=complex)
    if name == "s":
        return np.diag([1, 1j])
    if name == "sdg":
        return np.diag([1, -1j])
    if name in ("rx", "ry", "rz"):
        t = angle(gate["theta"])
        c, s = math.cos(t / 2), math.sin(t / 2)
        if name == "rx":
            return np.array([[c, -1j * s], [-1j * s, c]])
        if name == "ry":
            return np.array([[c, -s], [s, c]])
        return np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])
    if name == "u3":
        t, p, l = (angle(gate[key]) for key in ("theta", "phi", "lam"))
        c, s = math.cos(t / 2), math.sin(t / 2)
        return np.array([[c, -np.exp(1j * l) * s],
                         [np.exp(1j * p) * s, np.exp(1j * (p + l)) * c]])
    raise ValueError(f"unsupported gate: {name!r}")


def embedded_one_qubit(matrix, qubit, n):
    factors = [matrix if j == qubit else np.eye(2) for j in range(n)]
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor)
    return result


def cnot(control, target, n):
    size = 1 << n
    result = np.zeros((size, size), dtype=complex)
    c_mask, t_mask = 1 << (n - 1 - control), 1 << (n - 1 - target)
    for basis in range(size):
        output = basis ^ t_mask if basis & c_mask else basis
        result[output, basis] = 1
    return result


def right_shift(n):
    size = 1 << n
    result = np.zeros((size, size), dtype=complex)
    for basis in range(size):
        output = ((basis & 1) << (n - 1)) | (basis >> 1)
        result[output, basis] = 1
    return result


def total_spin_casimir(n):
    """Return S^2 up to an additive scalar (the sum of pairwise SWAPs)."""
    size = 1 << n
    result = np.zeros((size, size), dtype=complex)
    for first in range(n):
        for second in range(first + 1, n):
            first_mask = 1 << (n - 1 - first)
            second_mask = 1 << (n - 1 - second)
            for basis in range(size):
                if bool(basis & first_mask) == bool(basis & second_mask):
                    output = basis
                else:
                    output = basis ^ first_mask ^ second_mask
                result[output, basis] += 1
    return result


def evaluate(candidate, tolerance=1e-9):
    n = candidate["n"]
    if not isinstance(n, int) or isinstance(n, bool) or not 2 <= n <= 8:
        raise ValueError("n must be an integer from 2 to 8")
    gates = candidate["gates"]
    if not isinstance(gates, list):
        raise ValueError("gates must be a list")
    unitary = np.eye(1 << n, dtype=complex)
    cx_count = 0
    for position, gate in enumerate(gates):
        try:
            name = gate["gate"].lower()
            if name == "cx":
                control, target = gate["control"], gate["target"]
                if (not isinstance(control, int) or isinstance(control, bool)
                        or not isinstance(target, int) or isinstance(target, bool)
                        or control == target or not 0 <= control < n
                        or not 0 <= target < n):
                    raise ValueError("invalid CNOT control or target")
                matrix = cnot(control, target, n)
                cx_count += 1
            else:
                qubit = gate["qubit"]
                if not isinstance(qubit, int) or isinstance(qubit, bool) or not 0 <= qubit < n:
                    raise ValueError("invalid one-qubit target")
                matrix = embedded_one_qubit(one_qubit_matrix(gate), qubit, n)
            unitary = matrix @ unitary
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"gate {position}: {error}") from error
    diagonalized = unitary @ right_shift(n) @ unitary.conj().T
    diagonal = np.diag(diagonalized)
    off_diagonal = diagonalized - np.diag(diagonal)
    unitary_error = float(np.max(np.abs(unitary @ unitary.conj().T - np.eye(1 << n))))
    off_diagonal_error = float(np.max(np.abs(off_diagonal)))
    root_error = float(np.max(np.abs(diagonal ** n - 1)))
    spin_matrix = unitary @ total_spin_casimir(n) @ unitary.conj().T
    spin_error = float(np.max(np.abs(spin_matrix - np.diag(np.diag(spin_matrix)))))
    valid = max(unitary_error, off_diagonal_error, root_error) <= tolerance
    return {
        "n": n,
        "cnot_count": cx_count,
        "valid_diagonalizer": valid,
        "beats_18_cnot_reference": n == 4 and valid and cx_count < 18,
        "unitarity_error": unitary_error,
        "off_diagonal_error": off_diagonal_error,
        "eigenvalue_root_error": root_error,
        "total_spin_off_diagonal_error": spin_error,
        "diagonalizes_total_spin": spin_error <= tolerance,
        "diagonal_eigenvalues": [[float(z.real), float(z.imag)] for z in diagonal],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    result = evaluate(json.loads(args.candidate.read_text()), args.tolerance)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["valid_diagonalizer"] else 1)


if __name__ == "__main__":
    main()
