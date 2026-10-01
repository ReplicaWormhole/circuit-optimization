#!/usr/bin/env python3
"""One independently implemented exact candidate check in Q[z]/(z^16+1).

z=exp(i*pi/16); no NumPy, SymPy or repository circuit checker is imported.
Execution is allowed only from the matching reserved run workspace.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import fractions
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import sqlite3
import sys
import time

THREAD_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for key in THREAD_KEYS:
    os.environ[key] = "1"

HERE = Path(__file__).resolve().parent
FIXED = {
    "schema_version": 1,
    "purpose": "independent exact row-eigenvector and unitarity certificate of one complete 11-CX circuit",
    "method": "exact Fraction polynomial arithmetic in Q[z]/(z^16+1)",
    "seed": 0, "random_draws": 0, "thread_count": 1,
    "time_limit_seconds": 90, "candidate_filename": "candidate.json",
    "expected_n": 4, "expected_cx_count": 11, "board_id": 138,
    "root_order": ["1", "i", "-1", "-i"],
    "expected_multiplicities": [6, 3, 4, 3],
    "polynomial_coefficients": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "z_embedding": "exp(i*pi/16)",
    "evidence_scope": "Exact certificate for the saved literal ancilla-free gate list only; no optimality, total-spin, strong Schur-transform, or family exclusion claim.",
}
REQUIRED_KEYS = set(FIXED) | {"candidate_sha256", "script_sha256", "runtime_versions", "dependency_sha256"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def versions():
    return {"python": ".".join(map(str, sys.version_info[:3]))}


def dependencies():
    return {"python_executable": sha(sys.executable), "fractions": sha(fractions.__file__),
            "math_extension": sha(math.__file__)}


def write(path, payload):
    temporary = Path(path).with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


@dataclass(frozen=True)
class C:
    """Unique sixteen rational coefficients modulo z^16=-1."""
    a: tuple[Fraction, ...]

    def __add__(self, other):
        return C(tuple(x + y for x, y in zip(self.a, other.a)))

    def __neg__(self):
        return C(tuple(-x for x in self.a))

    def __sub__(self, other):
        return self + (-other)

    def scale(self, scalar):
        return C(tuple(x * scalar for x in self.a))

    def __mul__(self, other):
        out = [Fraction(0)] * 16
        for j, x in enumerate(self.a):
            if not x:
                continue
            for k, y in enumerate(other.a):
                if y:
                    degree = j + k
                    out[degree % 16] += x * y * (1 if degree < 16 else -1)
        return C(tuple(out))

    def conj(self):
        out = ZERO
        for exponent, coefficient in enumerate(self.a):
            if coefficient:
                out = out + zpower(-exponent).scale(coefficient)
        return out

    def serialize(self):
        return [[coefficient.numerator, coefficient.denominator] for coefficient in self.a]


ZERO = C((Fraction(0),) * 16)


def zpower(exponent):
    exponent %= 32
    out = [Fraction(0)] * 16
    out[exponent % 16] = Fraction(1 if exponent < 16 else -1)
    return C(tuple(out))


ONE = zpower(0)
I = zpower(8)


ANGLE = re.compile(r"^([+-]?)(?:(\d+)\*?)?pi(?:/(\d+))?$")


def pi_fraction(raw):
    # Float angles are rejected: certificates must preserve exact input.
    if raw == 0 and type(raw) is int:
        return Fraction(0)
    if not isinstance(raw, str):
        raise ValueError(f"angle must be a rational-pi string: {raw!r}")
    match = ANGLE.fullmatch(raw.replace(" ", ""))
    if not match:
        raise ValueError(f"unsupported exact angle: {raw!r}")
    sign, numerator, denominator = match.groups()
    return Fraction((-1 if sign == "-" else 1) * int(numerator or 1), int(denominator or 1))


def phase(raw):
    exponent = 16 * pi_fraction(raw)
    if exponent.denominator != 1:
        raise ValueError("phase is outside the conductor32 ring")
    return zpower(exponent.numerator)


def cosine_sine(raw):
    exponent = 8 * pi_fraction(raw)
    if exponent.denominator != 1:
        raise ValueError("half-angle is outside the conductor32 ring")
    positive, negative = zpower(exponent.numerator), zpower(-exponent.numerator)
    return (positive + negative).scale(Fraction(1, 2)), ((positive - negative) * (-I)).scale(Fraction(1, 2))


def local(gate):
    name = gate["gate"]
    if name == "h":
        r = (zpower(4) + zpower(-4)).scale(Fraction(1, 2))
        return ((r, r), (r, -r))
    if name == "x":
        return ((ZERO, ONE), (ONE, ZERO))
    if name in {"s", "sdg"}:
        return ((ONE, ZERO), (ZERO, I if name == "s" else -I))
    if name in {"rx", "ry", "rz"}:
        c, s = cosine_sine(gate["theta"])
        if name == "rx":
            return ((c, (-I) * s), ((-I) * s, c))
        if name == "ry":
            return ((c, -s), (s, c))
        exponent = 8 * pi_fraction(gate["theta"])
        return ((zpower(-exponent.numerator), ZERO), (ZERO, zpower(exponent.numerator)))
    if name == "u3":
        c, s = cosine_sine(gate["theta"])
        p, l = phase(gate["phi"]), phase(gate["lam"])
        return ((c, -(l * s)), (p * s, p * l * c))
    raise ValueError(f"unsupported elementary gate {name!r}")


def valid_qubit(q):
    return type(q) is int and 0 <= q < 4


def build(candidate):
    if candidate["n"] != 4 or type(candidate["n"]) is not int:
        raise ValueError("certificate is specialized to exactly four qubits")
    U = [[ONE if row == col else ZERO for col in range(16)] for row in range(16)]
    cx_count = 0
    for index, gate in enumerate(candidate["gates"]):
        if gate["gate"] == "cx":
            c, t = gate["control"], gate["target"]
            if not valid_qubit(c) or not valid_qubit(t) or c == t:
                raise ValueError(f"invalid CX at gate {index}")
            cmask, tmask = 1 << (3 - c), 1 << (3 - t)
            # A CX is its own inverse; left multiplication permutes output rows.
            U = [U[row ^ tmask] if row & cmask else U[row] for row in range(16)]
            cx_count += 1
        else:
            q = gate["qubit"]
            if not valid_qubit(q):
                raise ValueError(f"invalid target at gate {index}")
            M = local(gate)
            mask = 1 << (3 - q)
            for row0 in range(16):
                if row0 & mask:
                    continue
                row1 = row0 | mask
                old0, old1 = U[row0], U[row1]
                U[row0] = [M[0][0] * a + M[0][1] * b for a, b in zip(old0, old1)]
                U[row1] = [M[1][0] * a + M[1][1] * b for a, b in zip(old0, old1)]
    return U, cx_count


def right_shift_column(basis):
    # Independent literal bit-label construction of |x3 x0 x1 x2>.
    x = [(basis >> (3 - q)) & 1 for q in range(4)]
    y = [x[3], x[0], x[1], x[2]]
    return sum(bit << (3 - q) for q, bit in enumerate(y))


def certificate(U):
    permutation = [right_shift_column(basis) for basis in range(16)]
    if sorted(permutation) != list(range(16)):
        raise AssertionError("right shift is not a permutation")
    roots = [zpower(8 * label) for label in range(4)]
    labels, row_failures = [], []
    for row in range(16):
        matches = [label for label, eigenvalue in enumerate(roots)
                   if all(U[row][permutation[col]] == eigenvalue * U[row][col] for col in range(16))]
        if len(matches) != 1:
            row_failures.append({"row": row, "matching_roots": matches})
            labels.append(None)
        else:
            labels.append(matches[0])
    conjugates = [[value.conj() for value in row] for row in U]
    unitary_failures = []
    for row in range(16):
        for other in range(16):
            inner = ZERO
            for col in range(16):
                inner = inner + U[row][col] * conjugates[other][col]
            expected = ONE if row == other else ZERO
            if inner != expected:
                unitary_failures.append({"row": row, "other_row": other,
                                         "residual": (inner - expected).serialize()})
    return {"right_shift_column_permutation": permutation, "row_root_labels": labels,
            "root_order": FIXED["root_order"], "multiplicities": [labels.count(k) for k in range(4)],
            "row_eigenvector_failures": row_failures, "unitarity_failures": unitary_failures,
            "exact_row_eigenvector_identity": not row_failures,
            "exact_unitarity": not unitary_failures}


def guard():
    if not HERE.name.isdigit() or HERE.parent.name != "runs" or HERE.parent.parent.name != "experiments":
        raise RuntimeError("execute only from experiments/runs/<reserved-id>")
    root, run_id = HERE.parents[2], int(HERE.name)
    cfg = json.loads((HERE / "config.json").read_text())
    if set(cfg) != REQUIRED_KEYS:
        raise RuntimeError(f"config-key mismatch: {sorted(set(cfg) ^ REQUIRED_KEYS)}")
    for key, expected in FIXED.items():
        if type(cfg[key]) is not type(expected) or cfg[key] != expected:
            raise RuntimeError(f"frozen configuration mismatch for {key}")
    if cfg["script_sha256"] != sha(__file__):
        raise RuntimeError("script hash mismatch")
    candidate_path = HERE / cfg["candidate_filename"]
    if cfg["candidate_sha256"] != sha(candidate_path):
        raise RuntimeError("candidate hash mismatch")
    if cfg["runtime_versions"] != versions() or cfg["dependency_sha256"] != dependencies():
        raise RuntimeError("runtime version/dependency hash mismatch")
    with sqlite3.connect(f"file:{root / 'experiments.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute("SELECT status,code_sha256,config_json,method FROM attempts WHERE id=?", (run_id,)).fetchone()
    if row is None or row[0] != "running" or row[1] != cfg["script_sha256"] or json.loads(row[2]) != cfg or row[3] != cfg["method"]:
        raise RuntimeError("active ledger code/config/method mismatch")
    with sqlite3.connect(f"file:{root / 'collaboration/board.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute("SELECT status,owner,run_id FROM hypotheses WHERE id=?", (cfg["board_id"],)).fetchone()
    if row != ("running", "numerical_sol", run_id):
        raise RuntimeError("active board ownership/run link mismatch")
    if any(os.environ[key] != "1" for key in THREAD_KEYS):
        raise RuntimeError("thread environment mismatch")
    return run_id, cfg, candidate_path


def main():
    started = time.monotonic()
    run_id, cfg, candidate_path = guard()
    with (HERE / "execution_started.json").open("x") as marker:
        json.dump({"run_id": run_id, "config_sha256": sha(HERE / "config.json"),
                   "utc_start": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, marker)
    remaining = cfg["time_limit_seconds"] - (time.monotonic() - started)
    if remaining <= 0:
        raise TimeoutError("preflight exhausted budget")
    def stop(signum, frame):
        raise TimeoutError("90-second exact-check budget reached")
    signal.signal(signal.SIGALRM, stop)
    signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        # Ring sanity checks also run only after reservation, before circuit work.
        if zpower(16) != -ONE or zpower(32) != ONE or zpower(1) * zpower(1).conj() != ONE:
            raise AssertionError("cyclotomic quotient/conjugation sanity check failed")
        candidate = json.loads(candidate_path.read_text())
        U, cx_count = build(candidate)
        write(HERE / "unitary_exact.json", {"coefficient_order": "z^0 through z^15",
                                            "entries": [[v.serialize() for v in row] for row in U]})
        result = certificate(U)
        valid = result["exact_row_eigenvector_identity"] and result["exact_unitarity"] and cx_count == cfg["expected_cx_count"] and result["multiplicities"] == cfg["expected_multiplicities"]
        result.update(status="completed", run_id=run_id, exact_valid_diagonalizer=valid,
                      cx_count=cx_count, gate_count=len(candidate["gates"]),
                      candidate_sha256=cfg["candidate_sha256"], script_sha256=cfg["script_sha256"],
                      config_sha256=sha(HERE / "config.json"), unitary_sha256=sha(HERE / "unitary_exact.json"),
                      runtime_versions=versions(), dependency_sha256=dependencies(),
                      runtime_module_paths={"fractions": fractions.__file__, "math": math.__file__, "python": sys.executable},
                      thread_environment={key: os.environ[key] for key in THREAD_KEYS},
                      seconds=time.monotonic() - started, evidence_scope=cfg["evidence_scope"])
        write(HERE / "result.json", result)
        print(json.dumps({key: result[key] for key in ("run_id", "exact_valid_diagonalizer", "exact_unitarity", "cx_count", "multiplicities", "row_root_labels", "seconds")}, sort_keys=True))
        return 0 if valid else 1
    except BaseException as error:
        write(HERE / "failure.json", {"run_id": run_id, "error_type": type(error).__name__,
                                      "error": str(error), "seconds": time.monotonic() - started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    raise SystemExit(main())
