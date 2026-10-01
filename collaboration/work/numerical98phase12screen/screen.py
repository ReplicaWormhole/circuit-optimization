#!/usr/bin/env python3
"""Frozen finite screen of 16 literal phase-aware 12-CX circuits.

Preparation is read-only with respect to circuit evaluation: this program must
only be run from a root-reserved experiments/runs/<ID> workspace after review.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

# Set thread limits before NumPy loads a BLAS implementation.
for _thread_env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_env] = "1"

import numpy as np

HERE = Path(__file__).resolve().parent
CONFIG_PATH = HERE / "config.json"
RESULT_PATH = HERE / "result.json"
VARIANT_DIR = HERE / "variants"
SOURCE_FILES = {"check_circuit.py": "check_circuit.py", "exact_check.py": "exact_check.py"}
REQUIRED_KEYS = {
    "schema_version", "purpose", "script_sha256", "source_sha256",
    "variant_count", "expected_cx_count", "variant_rule", "time_limit_seconds",
    "thread_count", "numeric_tolerance", "exact_check_timeout_seconds",
    "evidence_scope",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cx(control: int, target: int) -> dict:
    return {"gate": "cx", "control": control, "target": target}


def one(name: str, qubit: int, **kwargs) -> dict:
    return {"gate": name, "qubit": qubit, **kwargs}


def h(q: int) -> dict:
    return one("h", q)


def x(q: int) -> dict:
    return one("x", q)


def ry(q: int, theta: str) -> dict:
    return one("ry", q, theta=theta)


def phase(q: int, theta: str) -> dict:
    # U3(0,0,theta) is diag(1, exp(i theta)); pi/4 is T, -pi/4 is Tdg.
    return one("u3", q, theta="0*pi", phi="0*pi", lam=theta)


def inverse_gate(gate: dict) -> dict:
    result = dict(gate)
    name = result["gate"]
    if name == "u3":
        result["lam"] = "-" + result["lam"] if not result["lam"].startswith("-") else result["lam"][1:]
    elif name == "ry":
        result["theta"] = "-" + result["theta"] if not result["theta"].startswith("-") else result["theta"][1:]
    elif name == "s":
        result["gate"] = "sdg"
    elif name == "sdg":
        result["gate"] = "s"
    elif name not in {"h", "x", "cx"}:
        raise ValueError(f"no inverse rule for {name}")
    return result


def rccx(control_a: int, control_b: int, target: int, adjoint: bool) -> list[dict]:
    # Frozen chronological literal from the task: H,T,CX_b,Tdg,CX_a,T,CX_b,Tdg,H.
    gates = [
        h(target), phase(target, "pi/4"), cx(control_b, target),
        phase(target, "-pi/4"), cx(control_a, target), phase(target, "pi/4"),
        cx(control_b, target), phase(target, "-pi/4"), h(target),
    ]
    if not adjoint:
        return gates
    return [inverse_gate(gate) for gate in reversed(gates)]


def make_variant(p_swap: bool, p_adjoint: bool,
                 cch_swap: bool, cch_adjoint: bool) -> dict:
    gates: list[dict] = []

    # Exact Bell decoder E_pair^dagger: CX then H on (0,2) and (1,3).
    gates += [cx(0, 2), h(0), cx(1, 3), h(1)]
    # Difference register R.
    gates += [cx(0, 1), cx(2, 3)]
    # Btilde=CZ_A, the 1-CX H-target-CX-H-target replacement.
    gates += [h(2), cx(0, 2), h(2)]

    # P: relative-phase Toffoli, choices (q3,q0) or swapped, forward/adjoint.
    p_controls = (0, 3) if p_swap else (3, 0)
    gates += rccx(p_controls[0], p_controls[1], 2, p_adjoint)

    # CH(q1 -> q0), exact one-CX implementation.
    gates += [ry(0, "-pi/4"), h(0), cx(1, 0), h(0), ry(0, "pi/4")]

    # Negative control q1=0 and positive control q3=1.  Chronologically apply
    # V^dagger=Ry(pi/4),S; RCCX; V=Sdg,Ry(-pi/4), with V on target q2.
    gates += [x(1), ry(2, "pi/4"), one("s", 2)]
    cch_controls = (3, 1) if cch_swap else (1, 3)
    gates += rccx(cch_controls[0], cch_controls[1], 2, cch_adjoint)
    gates += [one("sdg", 2), ry(2, "-pi/4"), x(1)]

    cx_count = sum(gate["gate"] == "cx" for gate in gates)
    if cx_count != 12:
        raise RuntimeError(f"variant has {cx_count} CX gates, expected 12")
    variant_id = f"p_swap{int(p_swap)}_p_adj{int(p_adjoint)}_cch_swap{int(cch_swap)}_cch_adj{int(cch_adjoint)}"
    return {
        "n": 4,
        "gates": gates,
        "metadata": {
            "variant_id": variant_id,
            "description": "Bell decoder + difference router + CZ_A + two literal RCCX selector variants",
            "evidence_scope": "one of exactly 16 predeclared 12-CX phase-aware variants; not an incumbent",
            "choices": {"p_control_swap": p_swap, "p_adjoint": p_adjoint,
                        "cch_control_swap": cch_swap, "cch_adjoint": cch_adjoint},
        },
    }


def make_variants() -> list[dict]:
    return [make_variant(*bits) for bits in itertools.product((False, True), repeat=4)]


def root_path() -> Path:
    if HERE.parent.name != "runs" or HERE.parent.parent.name != "experiments":
        raise RuntimeError("screen may execute only from experiments/runs/<reserved-id>/")
    root = HERE.parents[2]
    if not (root / "experiments.sqlite3").is_file():
        raise RuntimeError("experiment ledger not found")
    return root


def guard(root: Path, run_id: int) -> dict:
    config = json.loads(CONFIG_PATH.read_text())
    if set(config) != REQUIRED_KEYS:
        raise RuntimeError(f"config-key mismatch: {sorted(set(config) ^ REQUIRED_KEYS)}")
    if config["schema_version"] != 1 or config["purpose"] != "one deterministic finite screen of 16 literal phase-aware 12-CX variants":
        raise RuntimeError("configuration identity mismatch")
    if sha256(Path(__file__).resolve()) != config["script_sha256"]:
        raise RuntimeError("script hash mismatch")
    if {name: sha256(root / rel) for name, rel in SOURCE_FILES.items()} != config["source_sha256"]:
        raise RuntimeError("checker source hash mismatch")
    if config["variant_count"] != 16 or config["expected_cx_count"] != 12:
        raise RuntimeError("variant/cost configuration mismatch")
    if config["variant_rule"] != "cartesian product of P control swap/adjoint and CCH control swap/adjoint bits, False then True":
        raise RuntimeError("variant enumeration mismatch")
    if config["time_limit_seconds"] != 30 or config["thread_count"] != 1:
        raise RuntimeError("time/thread budget mismatch")
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        row = db.execute("SELECT status, code_sha256 FROM attempts WHERE id=?", (run_id,)).fetchone()
    if row is None or row[0] != "running" or row[1] != config["script_sha256"]:
        raise RuntimeError("run is not the matching active ledger reservation")
    return config


def evaluate_candidate(candidate: dict, tolerance: float) -> dict:
    # Run the repository checker for its matrix-level validity metrics, then
    # independently obtain the same full matrix through its documented helpers
    # to record normalized Frobenius off-diagonal loss as well as max error.
    import check_circuit

    checked = check_circuit.evaluate(candidate, tolerance=tolerance)
    unitary = np.eye(16, dtype=complex)
    cx_count = 0
    for gate in candidate["gates"]:
        if gate["gate"] == "cx":
            matrix = check_circuit.cnot(gate["control"], gate["target"], 4)
            cx_count += 1
        else:
            local = check_circuit.one_qubit_matrix(gate)
            matrix = check_circuit.embedded_one_qubit(local, gate["qubit"], 4)
        unitary = matrix @ unitary
    target = check_circuit.right_shift(4)
    conjugated = unitary @ target @ unitary.conj().T
    off = conjugated - np.diag(np.diag(conjugated))
    roots = (1 + 0j, 1j, -1 + 0j, -1j)
    labels = [min(range(4), key=lambda k: abs(complex(value[0], value[1]) - roots[k]))
              for value in checked["diagonal_eigenvalues"]]
    return {
        "variant_id": candidate["metadata"]["variant_id"],
        "candidate_sha256": hashlib.sha256((json.dumps(candidate, sort_keys=True) + "\n").encode()).hexdigest(),
        "cx_count": cx_count,
        "loss": float(np.vdot(off, off).real / 16.0),
        "maxoff": float(np.max(np.abs(off))),
        "unitarity_error": checked["unitarity_error"],
        "root_error": checked["eigenvalue_root_error"],
        "valid_diagonalizer": checked["valid_diagonalizer"],
        "eigenvalue_multiplicities": [labels.count(label) for label in range(4)],
        "candidate": candidate,
    }


def main() -> int:
    run_id = int(HERE.name) if HERE.name.isdigit() else -1
    root = root_path()
    config = guard(root, run_id)
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        if env[key] != "1":
            raise RuntimeError("single-thread environment guard failed")

    deadline = time.monotonic() + config["time_limit_seconds"] - 2.0
    VARIANT_DIR.mkdir(exist_ok=True)
    records = []
    started = time.monotonic()
    for index, candidate in enumerate(make_variants()):
        if time.monotonic() >= deadline:
            break
        path = VARIANT_DIR / f"{index:02d}_{candidate['metadata']['variant_id']}.json"
        path.write_text(json.dumps(candidate, indent=2, sort_keys=True) + "\n")
        record = evaluate_candidate(candidate, config["numeric_tolerance"])
        record["candidate_file"] = path.name
        records.append(record)
        (HERE / "partial_result.json").write_text(json.dumps({
            "run_id": run_id, "script_sha256": config["script_sha256"],
            "evaluated_count": len(records), "records": records,
        }, indent=2, sort_keys=True) + "\n")

    if len(records) != config["variant_count"]:
        raise TimeoutError(f"completed only {len(records)} of 16 frozen variants")
    best = min(records, key=lambda row: (row["loss"], row["variant_id"]))
    best_candidate_path = HERE / "best_candidate.json"
    best_candidate_path.write_text(json.dumps(best["candidate"], indent=2, sort_keys=True) + "\n")
    (HERE / "all_metrics.json").write_text(json.dumps([
        {key: value for key, value in row.items() if key != "candidate"} for row in records
    ], indent=2, sort_keys=True) + "\n")

    exact = None
    if best["valid_diagonalizer"]:
        command = [sys.executable, str(root / "exact_check.py"), str(best_candidate_path)]
        try:
            completed = subprocess.run(command, cwd=root, capture_output=True, text=True,
                                       env=env, timeout=config["exact_check_timeout_seconds"])
            exact = {"returncode": completed.returncode, "stdout": completed.stdout,
                     "stderr": completed.stderr, "timed_out": False}
        except subprocess.TimeoutExpired as error:
            exact = {"returncode": None,
                     "stdout": error.stdout.decode() if isinstance(error.stdout, bytes) else (error.stdout or ""),
                     "stderr": error.stderr.decode() if isinstance(error.stderr, bytes) else (error.stderr or ""),
                     "timed_out": True}

    result = {
        "run_id": run_id,
        "script_sha256": config["script_sha256"],
        "config_sha256": sha256(CONFIG_PATH),
        "source_sha256": config["source_sha256"],
        "evaluated_count": len(records),
        "variant_count": 16,
        "cx_count_each": 12,
        "seconds": time.monotonic() - started,
        "numpy_version": np.__version__,
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "best": {key: value for key, value in best.items() if key != "candidate"},
        "records": [{key: value for key, value in row.items() if key != "candidate"} for row in records],
        "exact_check": exact,
        "evidence_scope": config["evidence_scope"],
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
