"""Test output-eigenspace rotations against the fixed Fourier-boundary rank bound.

For the exact baseline U=S B P, rotate two output basis states with the same
V4 eigenvalue by W. Then U'=W U still diagonalizes V4 and its induced
boundary is B'=S^dagger W S B. This experiment tests whether B' can have
operator-Schmidt rank at most eight across 01|23, a necessary condition for
three crossing CNOTs. It does not synthesize W or B'.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import cnot, embedded_one_qubit, one_qubit_matrix
from exact_check import exact_eigenvalue_labels
from topology14_boundary_schmidt import matrix as boundary_matrix
from topology14_boundary_schmidt import schmidt_singular_values


ROOT = Path(__file__).resolve().parent


def circuit_matrix(gates):
    out = np.eye(16, dtype=complex)
    for gate in gates:
        block = (cnot(gate["control"], gate["target"], 4)
                 if gate["gate"] == "cx"
                 else embedded_one_qubit(one_qubit_matrix(gate), gate["qubit"], 4))
        out = block @ out
    return out


def two_level(first, second, theta, axis):
    out = np.eye(16, dtype=complex)
    cosine, sine = math.cos(theta), math.sin(theta)
    if axis == "y":
        offdiag = (-sine, sine)
    elif axis == "x":
        offdiag = (-1j * sine, -1j * sine)
    else:
        raise ValueError(axis)
    out[first, first] = out[second, second] = cosine
    out[first, second], out[second, first] = offdiag
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    labels = exact_eigenvalue_labels(baseline)["output_labels"]
    suffix = circuit_matrix(baseline["gates"][38:])
    boundary = boundary_matrix(0.0)
    initial_singular = schmidt_singular_values(boundary, (0, 1))
    assert np.min(initial_singular) > 1e-8

    records = []
    for first in range(16):
        for second in range(first + 1, 16):
            if labels[first] != labels[second]:
                continue
            for axis in ("x", "y"):
                for numerator in (1, 2, 3, 4):
                    theta = numerator * math.pi / 8
                    w = two_level(first, second, theta, axis)
                    changed = suffix.conj().T @ w @ suffix @ boundary
                    singular = schmidt_singular_values(changed, (0, 1))
                    records.append({"output_pair": [first, second],
                                    "label": labels[first], "axis": axis,
                                    "angle_over_pi": f"{numerator}/8",
                                    "rank_1e-9": int(np.count_nonzero(singular > 1e-9)),
                                    "smallest_singular": float(singular[-1])})
    minimum = min(record["smallest_singular"] for record in records)
    rank_hist = {str(rank): sum(r["rank_1e-9"] == rank for r in records)
                 for rank in sorted({r["rank_1e-9"] for r in records})}
    result = {"trial_count": len(records), "rank_histogram": rank_hist,
              "smallest_singular_over_trials": minimum,
              "best": sorted(records, key=lambda r: r["smallest_singular"])[:8],
              "necessary_condition": "rank <= 8 for a three-CNOT boundary across 01|23",
              "scope": "one two-level eigenspace SU2 rotation; 2 axes and four exact pi/8 angles"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
