"""Independent statevector check of saved distance-three fit candidates."""

import json
import math
from pathlib import Path

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "search13_fulltail_rank8_distance3_fit_result.json"
OUT = ROOT / "search13_fulltail_rank8_distance3_validation.json"


def angle(value):
    return float(sp.N(sp.sympify(value, locals={"pi": sp.pi})))


def local_matrix(gate):
    name = gate["gate"]
    if name == "h":
        return np.array([[1, 1], [1, -1]], complex) / math.sqrt(2)
    if name == "x":
        return np.array([[0, 1], [1, 0]], complex)
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
        return np.diag([complex(math.cos(t / 2), -math.sin(t / 2)),
                        complex(math.cos(t / 2), math.sin(t / 2))])
    if name == "u3":
        t, p, l = (angle(gate[k]) for k in ("theta", "phi", "lam"))
        c, s = math.cos(t / 2), math.sin(t / 2)
        return np.array([[c, -np.exp(1j * l) * s],
                         [np.exp(1j * p) * s, np.exp(1j * (p + l)) * c]])
    raise ValueError(name)


def circuit_matrix(gates):
    states = np.eye(16, dtype=complex)
    for gate in gates:
        if gate["gate"] == "cx":
            cmask = 1 << (3 - gate["control"])
            tmask = 1 << (3 - gate["target"])
            for row in range(16):
                if row & cmask and not row & tmask:
                    states[[row, row | tmask], :] = states[[row | tmask, row], :]
            continue
        mask = 1 << (3 - gate["qubit"])
        local = local_matrix(gate)
        for row in range(16):
            if not row & mask:
                pair = states[[row, row | mask], :].copy()
                states[[row, row | mask], :] = local @ pair
    return states


def main():
    assert not OUT.exists()
    fit = json.loads(RESULT.read_text())
    shift = np.zeros((16, 16), complex)
    for basis in range(16):
        shifted = ((basis & 1) << 3) | (basis >> 1)
        shift[shifted, basis] = 1
    rows = []
    for case in fit["rows"]:
        for stage in ("process", "direct", "warm_direct"):
            info = case[stage]
            candidate = json.loads((ROOT / info["candidate"]).read_text())
            gates = candidate["gates"]
            matrix = circuit_matrix(gates)
            diagonalized = matrix @ shift @ matrix.conj().T
            off = diagonalized - np.diag(np.diag(diagonalized))
            error = float(np.max(np.abs(off)))
            loss = float(np.sum(np.abs(off) ** 2).real / 16)
            row = {"case": case["case"], "stage": stage,
                   "candidate": info["candidate"],
                   "cnot_count": sum(g["gate"] == "cx" for g in gates),
                   "off_diagonal_error": error,
                   "label_free_loss": loss,
                   "checker_error_delta": abs(error - info["off_diagonal_error"]),
                   "unitarity_error": float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(16))))}
            if stage != "process":
                row["optimizer_loss_delta"] = abs(loss - info["loss"])
            rows.append(row)
    assert len(rows) == 18
    assert all(r["cnot_count"] == 13 and r["off_diagonal_error"] > 1e-3
               and r["checker_error_delta"] < 1e-11 and r["unitarity_error"] < 1e-11
               and r.get("optimizer_loss_delta", 0) < 1e-11 for r in rows)
    OUT.write_text(json.dumps({"source": RESULT.name, "rows": rows,
                               "method": "independent basis-state amplitude simulation; no checker or optimizer imports"},
                              indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidate_count": len(rows),
                      "min_offdiag": min(r["off_diagonal_error"] for r in rows),
                      "min_direct_loss": min(r["label_free_loss"] for r in rows),
                      "max_checker_delta": max(r["checker_error_delta"] for r in rows),
                      "max_optimizer_delta": max(r.get("optimizer_loss_delta", 0) for r in rows)}))


if __name__ == "__main__":
    main()
