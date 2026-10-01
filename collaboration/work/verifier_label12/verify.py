"""Independent gate-list matrix check; no shared simulation imports."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


def cycle():
    v = np.zeros((16, 16), complex)
    for x in range(16):
        v[(x >> 1) | ((x & 1) << 3), x] = 1
    return v


def matrix(data):
    if data["n"] != 4:
        raise ValueError("Four qubits required")
    u = np.eye(16, dtype=complex)
    count = 0
    for g in data["gates"]:
        if g["gate"] == "cx":
            c, t = g["control"], g["target"]
            if c == t or c not in range(4) or t not in range(4):
                raise ValueError("Invalid CX wires")
            op = np.zeros((16, 16), complex)
            for x in range(16):
                y = x ^ (1 << (3-t)) if (x >> (3-c)) & 1 else x
                op[y, x] = 1
            count += 1
        elif g["gate"] == "u3":
            q = g["qubit"]
            if q not in range(4):
                raise ValueError("Invalid local wire")
            a, b, d = [g[k] for k in ("theta", "phi", "lam")]
            cs, sn = np.cos(a/2), np.sin(a/2)
            local = np.array([[cs, -np.exp(1j*d)*sn],
                              [np.exp(1j*b)*sn, np.exp(1j*(b+d))*cs]])
            op = np.array([[1]], complex)
            for wire in range(4):
                op = np.kron(op, local if wire == q else np.eye(2))
        else:
            raise ValueError(f"Unsupported gate {g['gate']}")
        u = op @ u
    return u, count


def check(path):
    p = Path(path)
    u, count = matrix(json.loads(p.read_text()))
    v = cycle()
    w = u @ v @ u.conj().T
    diag = np.diag(w)
    off = w - np.diag(diag)
    roots = np.array([1, 1j, -1, -1j])
    labels = np.argmin(abs(diag[:, None] - roots), axis=1)
    d = np.diag(roots[labels])
    unitarity = float(np.max(abs(u.conj().T @ u - np.eye(16))))
    maxoff = float(np.max(abs(off)))
    rooterror = float(np.max(abs(diag-roots[labels])))
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "cx_count": count, "unitarity_error": unitarity,
            "max_offdiagonal": maxoff, "offdiagonal_loss": float(np.sum(abs(off)**2)/16),
            "nearest_root_max_error": rooterror, "nearest_labels": labels.tolist(),
            "nearest_label_counts": np.bincount(labels, minlength=4).tolist(),
            "nearest_label_intertwining_max_error": float(np.max(abs(u@v-d@u))),
            "numerically_valid": bool(unitarity < 1e-10 and maxoff < 1e-10 and rooterror < 1e-10),
            "exact_certificate": False}


if __name__ == "__main__":
    report = {"convention": "q0 MSB; chronological left multiplication; right-cycle UV=DU",
              "cycle_traces": [int(np.trace(np.linalg.matrix_power(cycle(), k)).real) for k in range(4)],
              "required_root_counts": [6, 3, 4, 3], "checks": [check(p) for p in sys.argv[2:]]}
    Path(sys.argv[1]).write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
