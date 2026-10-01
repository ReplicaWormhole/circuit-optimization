"""Scan exact mid-layer Z phases after deleting exact15's middle CZ.

Keep the existing exact rational-pi two-qubit pair blocks and prefix.
Delete chronological H1-CX(2,1)-H1 at gates 49:52, then insert
independent Rz(k*pi/denominator) on each wire between the two pair layers.
This finite restricted family tests whether simple mid-layer phases alone
replace the CZ while choosing a different V4 eigenbasis.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from check_circuit import evaluate, right_shift
from topology14_twolevel_boundary_rank import circuit_matrix


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"


def phase_vector(values, denominator):
    angles = np.array(values, dtype=float) * math.pi / denominator
    out = np.empty(16, dtype=complex)
    for x in range(16):
        z = np.array([1 - 2 * ((x >> (3 - q)) & 1) for q in range(4)])
        out[x] = np.exp(-0.5j * np.dot(angles, z))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--denominator", type=int, required=True)
    args = parser.parse_args()
    denominator = args.denominator
    assert denominator >= 1
    source = json.loads(SOURCE.read_text())
    gates = source["gates"]
    assert gates[49:52] == [{"gate": "h", "qubit": 1},
                            {"gate": "cx", "control": 2, "target": 1},
                            {"gate": "h", "qubit": 1}]
    prefix_first = circuit_matrix(gates[:49])
    final = circuit_matrix(gates[52:])
    target = prefix_first @ right_shift(4) @ prefix_first.conj().T
    best = None
    counts = {}
    for values in itertools.product(range(2 * denominator), repeat=4):
        phases = phase_vector(values, denominator)
        changed = final @ ((phases[:, None] * target) * phases.conj()[None, :]) @ final.conj().T
        offdiag = changed - np.diag(np.diag(changed))
        error = float(np.max(np.abs(offdiag)))
        bucket = round(error, 6)
        counts[bucket] = counts.get(bucket, 0) + 1
        if best is None or error < best["off_diagonal_error"]:
            best = {"phase_numerators": values, "off_diagonal_error": error}
    inserted = [{"gate": "rz", "qubit": q,
                 "theta": f"{value}*pi/{denominator}"}
                for q, value in enumerate(best["phase_numerators"])]
    candidate = {"n": 4, "gates": gates[:49] + inserted + gates[52:]}
    name = f"exact14_mid_z_grid_{denominator}"
    best_path = ROOT / f"{name}_best.json"
    best_path.write_text(json.dumps(candidate, indent=2) + "\n")
    checked = evaluate(candidate)
    result = {"denominator": denominator,
              "trials": (2 * denominator) ** 4,
              "best": best, "checked_error": checked["off_diagonal_error"],
              "valid": checked["valid_diagonalizer"],
              "cnot_count": checked["cnot_count"],
              "near_zero_count": sum(n for value, n in counts.items()
                                     if value < 1e-6)}
    (ROOT / f"{name}_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
