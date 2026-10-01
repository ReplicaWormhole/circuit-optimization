#!/usr/bin/env python3
"""Reserved numerical and exact check of a 12-CX parity-difference baseline."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import time

for _thread_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_key] = "1"

import numpy as np
import sympy


WORK = Path(__file__).resolve().parent
CONFIG_PATH = WORK / "config.json"
RESULT_PATH = WORK / "result.json"
CANDIDATE_PATH = WORK / "candidate.json"
SOURCE_FILES = {
    "native_gate_check.py": "native_gate_check.py",
    "check_circuit.py": "check_circuit.py",
    "exact_check.py": "exact_check.py",
}
REQUIRED_CONFIG_KEYS = {
    "schema_version", "purpose", "script_sha256", "source_sha256",
    "runtime_versions", "expected_cnot_count", "expected_native_entanglers",
    "numerical_timeout_seconds", "exact_timeout_seconds", "thread_count",
    "numerical_command", "exact_command", "candidate_filename",
    "evidence_scope", "candidate_sha256", "seed", "random_draws", "external_timeout_seconds", "expected_exact_labels",
}


def cx(control: int, target: int) -> dict:
    return {"gate": "cx", "control": control, "target": target}


def h(qubit: int) -> dict:
    return {"gate": "h", "qubit": qubit}


def x(qubit: int) -> dict:
    return {"gate": "x", "qubit": qubit}


def ry(qubit: int, theta: str) -> dict:
    return {"gate": "ry", "qubit": qubit, "theta": theta}


def phase(qubit: int, theta: str) -> dict:
    # u3(0, 0, theta) is exactly diag(1, exp(i theta)).
    return {"gate": "u3", "qubit": qubit, "theta": "0*pi",
            "phi": "0*pi", "lam": theta}


def rccx(control_a: int, control_b: int, target: int) -> list[dict]:
    """Literal relative-phase Toffoli, chronological I/I/Z/Y blocks."""
    return [h(target), phase(target, "pi/4"), cx(control_b, target),
            phase(target, "-pi/4"), cx(control_a, target), phase(target, "pi/4"),
            cx(control_b, target), phase(target, "-pi/4"), h(target)]


def rz(qubit: int, theta: str) -> dict:
    return {"gate": "rz", "qubit": qubit, "theta": theta}


def rx(qubit: int, theta: str) -> dict:
    return {"gate": "rx", "qubit": qubit, "theta": theta}


def sg(qubit: int) -> dict:
    return {"gate": "s", "qubit": qubit}


def sdg(qubit: int) -> dict:
    return {"gate": "sdg", "qubit": qubit}


def make_gates() -> list[dict]:
    gates: list[dict] = []

    # E_pair^dagger = H_c CX(c->t), chronologically CX then H, on (0,2),(1,3).
    gates += [cx(0, 2), h(0), cx(1, 3), h(1)]

    # Difference register d=a xor b: (a,b)->(a,d).
    gates += [cx(0, 1), cx(2, 3)]

    # Exact four-CX merge of B=CSdag and corrected P controls(x=q0,v=q3).
    # Target reflections N1=(-X+Y)/sqrt2,N2=N4=(-X+Z)/sqrt2,N3=-X.
    gates += [rz(2, "-3*pi/4"), cx(0, 2), rz(2, "3*pi/4")]
    gates += [ry(2, "3*pi/4"), cx(3, 2), ry(2, "-3*pi/4")]
    gates += [rz(2, "-pi"), cx(0, 2), rz(2, "pi")]  # controlled(-X).
    gates += [ry(2, "3*pi/4"), cx(3, 2), ry(2, "-3*pi/4")]
    gates += [phase(0, "pi/4")]  # restore relative control-sector phase.

    # Exact CH; active Ry(+pi/4) Z Ry(-pi/4)=H.
    gates += [ry(0, "-pi/4"), h(0), cx(1, 0), h(0), ry(0, "pi/4")]

    # Negative control c=1-u, positive v; Vprime maps Y->H and Z->Y.
    # Vprime=Ry(-3pi/4) Rx(-pi/2), apply its dagger before RCCX.
    gates += [x(1), ry(2, "3*pi/4"), rx(2, "pi/2")]
    gates += rccx(1, 3, 2)
    gates += [rx(2, "-pi/2"), ry(2, "-3*pi/4"), x(1)]

    return gates


def make_candidate() -> dict:
    return {
        "n": 4,
        "gates": make_gates(),
        "metadata": {
            "description": "Bell decoder, parity-difference router, common-phase repaired relative-phase controlled block diagonalizer",
            "evidence_scope": "12-CX candidate pending independent gate-list acceptance",
            "chronological_blocks": ["E_pair_dagger", "R", "merged_CSdag_corrected_RCCX_P_4CX", "controlled_H0", "Vprime_RCCX_CCH"],
        },
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_root() -> Path:
    for parent in WORK.parents:
        if (parent / "experiments.sqlite3").is_file() and (parent / "exact_check.py").is_file():
            return parent
    raise RuntimeError("could not identify repository root")


def guarded_inputs(root: Path, run_id: int) -> dict:
    config = json.loads(CONFIG_PATH.read_text())
    if set(config) != REQUIRED_CONFIG_KEYS:
        raise RuntimeError(f"config key mismatch: {sorted(set(config) ^ REQUIRED_CONFIG_KEYS)}")
    if config["schema_version"] != 1 or config["purpose"] != "one bounded numerical and exact check of an explicit 12-CX parity-difference baseline":
        raise RuntimeError("configuration identity mismatch")
    if sha256(Path(__file__).resolve()) != config["script_sha256"]:
        raise RuntimeError("frozen script hash mismatch")
    if config["runtime_versions"] != {
        "python": ".".join(map(str, sys.version_info[:3])),
        "numpy": np.__version__,
        "sympy": sympy.__version__,
    }:
        raise RuntimeError("runtime version mismatch")
    actual_sources = {name: sha256(root / rel) for name, rel in SOURCE_FILES.items()}
    if actual_sources != config["source_sha256"]:
        raise RuntimeError("checker source hash mismatch")
    if config["candidate_filename"] != CANDIDATE_PATH.name:
        raise RuntimeError("candidate filename mismatch")
    if config["seed"] != 0 or config["random_draws"] != 0:
        raise RuntimeError("deterministic seed/no-randomness guard failed")
    if config["thread_count"] != 1:
        raise RuntimeError("thread-count guard failed")
    if config["numerical_command"] != ["native_gate_check.py", "candidate.json", "--tolerance", "1e-9"]:
        raise RuntimeError("numerical command mismatch")
    if config["exact_command"] != ["exact_check.py", "candidate.json"]:
        raise RuntimeError("exact command mismatch")

    with sqlite3.connect(root / "experiments.sqlite3") as db:
        row = db.execute("SELECT status, code_sha256, config_json FROM attempts WHERE id=?", (run_id,)).fetchone()
    if row is None or row[0] != "running" or row[1] != config["script_sha256"]:
        raise RuntimeError("run is not the matching active ledger reservation")
    reserved = json.loads(row[2])
    if reserved.get("frozen_config") != config or reserved.get("frozen_config_sha256") != sha256(CONFIG_PATH):
        raise RuntimeError("ledger frozen config binding mismatch")
    if sha256(CANDIDATE_PATH) != config["candidate_sha256"]:
        raise RuntimeError("frozen candidate hash mismatch")
    if json.loads(CANDIDATE_PATH.read_text()) != make_candidate():
        raise RuntimeError("frozen candidate literal word mismatch")
    if config["numerical_timeout_seconds"] != 30 or config["exact_timeout_seconds"] != 60 or config["external_timeout_seconds"] != 100:
        raise RuntimeError("frozen timeout guard failed")
    if config["expected_exact_labels"] != [0,0,0,2,0,1,3,0,2,1,0,3,2,3,1,2]:
        raise RuntimeError("independent predicted-label guard mismatch")
    return config


def run_checker(command: list[str], timeout: int, root: Path, candidate: Path,
                env: dict[str, str]) -> dict:
    argv = [sys.executable, str(root / command[0]), str(candidate), *command[2:]]
    started = time.monotonic()
    try:
        completed = subprocess.run(argv, cwd=root, capture_output=True, text=True,
                                   env=env, timeout=timeout)
        return {"returncode": completed.returncode,
                "seconds": time.monotonic() - started,
                "timed_out": False, "stdout": completed.stdout,
                "stderr": completed.stderr}
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout.decode() if isinstance(error.stdout, bytes) else (error.stdout or "")
        stderr = error.stderr.decode() if isinstance(error.stderr, bytes) else (error.stderr or "")
        return {"returncode": None, "seconds": time.monotonic() - started,
                "timed_out": True, "stdout": stdout, "stderr": stderr}


def main() -> int:
    if not WORK.name.isdigit() or not WORK.parent.name == "runs" or not WORK.parent.parent.name == "experiments":
        raise RuntimeError("execution is permitted only from experiments/runs/<reserved-id>/")
    run_id = int(WORK.name)
    root = find_root()
    config = guarded_inputs(root, run_id)

    with (WORK / "execution_started.json").open("x") as handle:
        json.dump({"run_id": run_id, "script_sha256": config["script_sha256"],
                   "candidate_sha256": config["candidate_sha256"]}, handle)
    candidate = json.loads(CANDIDATE_PATH.read_text())
    cx_count = sum(gate["gate"] == "cx" for gate in candidate["gates"])
    if cx_count != config["expected_cnot_count"] or cx_count != 12:
        raise RuntimeError(f"unexpected CNOT count: {cx_count}")
    if any(gate["gate"] not in {"cx", "h", "x", "ry", "rx", "rz", "s", "sdg", "u3"}
           for gate in candidate["gates"]):
        raise RuntimeError("candidate contains a gate outside the frozen primitive list")
    if config["expected_native_entanglers"] != 0:
        raise RuntimeError("native-entangler count mismatch")
    candidate_hash = sha256(CANDIDATE_PATH)

    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    numerical = run_checker(config["numerical_command"], config["numerical_timeout_seconds"],
                            root, CANDIDATE_PATH, env)
    exact = run_checker(config["exact_command"], config["exact_timeout_seconds"],
                        root, CANDIDATE_PATH, env)

    (WORK / "numerical.stdout.json").write_text(numerical["stdout"])
    (WORK / "numerical.stderr.txt").write_text(numerical["stderr"])
    (WORK / "exact.stdout.json").write_text(exact["stdout"])
    (WORK / "exact.stderr.txt").write_text(exact["stderr"])
    try:
        numerical_pass = json.loads(numerical["stdout"])["valid_diagonalizer"] is True
        exact_info = json.loads(exact["stdout"])
        exact_pass = exact_info["exact_diagonalizer"] is True and exact_info["output_labels"] == config["expected_exact_labels"] and exact_info["cnot_count"] == 12
    except (ValueError, KeyError, TypeError):
        numerical_pass = exact_pass = False
    result = {
        "run_id": run_id,
        "script_sha256": config["script_sha256"],
        "config_sha256": sha256(CONFIG_PATH),
        "candidate_sha256": candidate_hash,
        "purpose": config["purpose"],
        "cnot_count": cx_count,
        "native_entanglers": config["expected_native_entanglers"],
        "expected_exact_labels": config["expected_exact_labels"],
        "numerical_returncode": numerical["returncode"],
        "numerical_seconds": numerical["seconds"],
        "numerical_timed_out": numerical["timed_out"],
        "exact_returncode": exact["returncode"],
        "exact_seconds": exact["seconds"],
        "exact_timed_out": exact["timed_out"],
        "status": "pass" if numerical["returncode"] == 0 and exact["returncode"] == 0
        and numerical_pass and exact_pass else "failed",
        "evidence_scope": config["evidence_scope"],
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
