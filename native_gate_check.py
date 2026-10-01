"""Numerical native-gate checker; MSB-first wires, chronological gates.

Native families: xx_yy(a,b)=exp(i*a*XX+i*b*YY), iswap=xx_yy(pi/4,pi/4),
sqrt_iswap=xx_yy(pi/8,pi/8), cp(theta)=diag(1,1,1,exp(i*theta)).
Angles are radians or the existing 'pi/8' notation. Costs are constructive
CNOT upper bounds with free arbitrary one-qubit gates, never optimality claims.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import angle, cnot, embedded_one_qubit, one_qubit_matrix, right_shift

NATIVE = {"xx_yy", "iswap", "sqrt_iswap", "cp"}
LOCAL = {"h", "x", "s", "sdg", "rx", "ry", "rz", "u3"}


def finite_angle(value):
    try:
        result = angle(value)
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as error:
        raise ValueError(f"invalid finite angle: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError("angle must be finite")
    return result


def wire(value, n):
    if not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < n:
        raise ValueError(f"invalid wire: {value!r}")
    return value


def native_matrix(gate):
    """Two-wire matrix in |00>,|01>,|10>,|11> order."""
    name = gate["gate"].lower()
    if name == "cp":
        return np.diag([1, 1, 1, np.exp(1j * finite_angle(gate["theta"]))])
    if name == "xx_yy":
        a, b = finite_angle(gate["a"]), finite_angle(gate["b"])
    elif name in {"iswap", "sqrt_iswap"}:
        a = b = math.pi / (4 if name == "iswap" else 8)
    else:
        raise ValueError(f"unsupported native gate: {name!r}")
    # XX and YY commute. Even/odd parity blocks rotate by a-b/a+b.
    even, odd = a - b, a + b
    return np.array([[math.cos(even), 0, 0, 1j * math.sin(even)],
                     [0, math.cos(odd), 1j * math.sin(odd), 0],
                     [0, 1j * math.sin(odd), math.cos(odd), 0],
                     [1j * math.sin(even), 0, 0, math.cos(even)]], complex)


def embedded_two_qubit(matrix, first, second, n):
    """Embed ordered local wires (first,second), including reversed wire order."""
    first, second = wire(first, n), wire(second, n)
    if first == second:
        raise ValueError("two-qubit wires must be distinct")
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError("two-qubit matrix must be finite and 4 by 4")
    result = np.zeros((1 << n, 1 << n), complex)
    masks = (1 << (n - 1 - first), 1 << (n - 1 - second))
    for col in range(1 << n):
        local_col = 2 * bool(col & masks[0]) + bool(col & masks[1])
        remaining = col & ~(masks[0] | masks[1])
        for local_row in range(4):
            row = remaining | (masks[0] if local_row & 2 else 0) | (masks[1] if local_row & 1 else 0)
            result[row, col] = matrix[local_row, local_col]
    return result


def circuit_matrix(candidate):
    n = candidate["n"]
    if not isinstance(n, int) or isinstance(n, bool) or not 2 <= n <= 8:
        raise ValueError("n must be an integer from 2 to 8")
    if not isinstance(candidate["gates"], list):
        raise ValueError("gates must be a list")
    result = np.eye(1 << n, dtype=complex)
    for position, gate in enumerate(candidate["gates"]):
        try:
            name = gate["gate"].lower()
            if name in NATIVE:
                pair = gate["qubits"]
                if not isinstance(pair, list) or len(pair) != 2:
                    raise ValueError("qubits must be a two-wire list")
                matrix = embedded_two_qubit(native_matrix(gate), *pair, n)
            elif name == "cx":
                control, target = wire(gate["control"], n), wire(gate["target"], n)
                if control == target:
                    raise ValueError("CNOT wires must be distinct")
                matrix = cnot(control, target, n)
            elif name in LOCAL:
                for key in ({"theta", "phi", "lam"} if name == "u3" else {"theta"} if name in {"rx", "ry", "rz"} else set()):
                    finite_angle(gate[key])
                matrix = embedded_one_qubit(one_qubit_matrix(gate), wire(gate["qubit"], n), n)
            else:
                raise ValueError(f"unsupported gate: {name!r}")
            if not np.isfinite(matrix).all() or np.max(abs(matrix.conj().T @ matrix - np.eye(len(matrix)))) > 1e-12:
                raise ValueError("gate matrix is not finite and unitary")
            result = matrix @ result
        except (KeyError, TypeError, ValueError, AttributeError, OverflowError) as error:
            raise ValueError(f"gate {position}: {error}") from error
    return result


def phase_error(first, second):
    """Max matrix error after the least-squares global phase alignment."""
    overlap = np.vdot(first, second)
    if abs(overlap) == 0:
        return float(np.max(abs(first - second)))
    return float(np.max(abs(first - second / (overlap / abs(overlap)))))


def costs(candidate):
    # Validate before reporting costs, so malformed inputs cannot get metrics.
    circuit_matrix(candidate)
    levels = [0] * candidate["n"]
    count = native_count = cx_count = 0
    compiled_upper_bound = 0
    for gate in candidate["gates"]:
        name = gate["gate"].lower()
        if name in NATIVE or name == "cx":
            pair = gate["qubits"] if name in NATIVE else [gate["control"], gate["target"]]
            layer = max(levels[q] for q in pair) + 1
            for q in pair:
                levels[q] = layer
            count += 1
            native_count += name in NATIVE
            cx_count += name == "cx"
            compiled_upper_bound += 1 if name == "cx" else 2
    return {"two_qubit_count": count, "native_family_gate_count": native_count,
            "cnot_count": cx_count, "two_qubit_depth": max(levels),
            "cnot_compilation_upper_bound": compiled_upper_bound,
            "minimum_cnot_cost": None}


def evaluate(candidate, tolerance=1e-9, reference=None):
    if not isinstance(tolerance, (float, int)) or isinstance(tolerance, bool) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    unitary = circuit_matrix(candidate)
    transformed = unitary @ right_shift(candidate["n"]) @ unitary.conj().T
    diagonal = np.diag(transformed)
    unitary_error = float(np.max(abs(unitary.conj().T @ unitary - np.eye(len(unitary)))))
    off_diagonal_error = float(np.max(abs(transformed - np.diag(diagonal))))
    root_error = float(np.max(abs(diagonal ** candidate["n"] - 1)))
    result = {**costs(candidate), "n": candidate["n"], "tolerance": tolerance,
              "unitarity_error": unitary_error, "off_diagonal_error": off_diagonal_error,
              "eigenvalue_root_error": root_error,
              "valid_diagonalizer": max(unitary_error, off_diagonal_error, root_error) <= tolerance,
              "proof_status": "numerical matrix evidence; this checker supplies no exact certificate"}
    if reference is not None:
        if candidate["n"] != reference["n"]:
            raise ValueError("reference n must match candidate n")
        result["reference_matrix_error_up_to_phase"] = phase_error(unitary, circuit_matrix(reference))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text()) if args.reference else None
    result = evaluate(json.loads(args.candidate.read_text()), args.tolerance, reference)
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result["valid_diagonalizer"] else 1)


if __name__ == "__main__":
    main()
