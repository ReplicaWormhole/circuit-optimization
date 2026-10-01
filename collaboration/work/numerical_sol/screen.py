#!/usr/bin/env python3
"""Prepare or once execute the frozen sixteen literal 12-CX words.

--prepare serializes gates and hashes only; it does not construct matrices.
Execution requires a matching active ledger row and board link in a run folder.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import signal
import sqlite3
import sys
import time

THREAD_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for key in THREAD_KEYS:
    os.environ[key] = "1"
import numpy as np

HERE = Path(__file__).resolve().parent
PURPOSE = "one deterministic finite NumPy screen of sixteen literal phase-aware 12-CX words"
RULE = "lexicographic False/True product: P control swap, P adjoint, CCH control swap, CCH adjoint"
SCOPE = "Finite literal-variant numerical evidence only; no optimizer, exact certificate, family exclusion, lower bound, or incumbent promotion. Adjoint choices may be unitary redundant."
FIXED = {
    "schema_version": 1, "purpose": PURPOSE, "seed": 0,
    "method": "complete NumPy unitary diagnostic", "variant_count": 16,
    "expected_cx_count": 12, "variant_rule": RULE, "time_limit_seconds": 30,
    "thread_count": 1, "numeric_tolerance": 1e-9,
    "root_order": ["1", "-1", "i", "-i"],
    "expected_multiplicities": [6, 4, 3, 3], "board_id": 120,
    "evidence_scope": SCOPE,
}
HASH_KEYS = {"script_sha256", "source_sha256", "dependency_sha256", "runtime_versions", "variants_sha256"}
REQUIRED = set(FIXED) | HASH_KEYS


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def cx(c, t):
    return {"gate": "cx", "control": c, "target": t}


def one(gate, q, **angles):
    return {"gate": gate, "qubit": q, **angles}


def phase(q, sign):
    return one("u3", q, theta="0*pi", phi="0*pi", lam="pi/4" if sign == 1 else "-pi/4")


def inverse(g):
    g = dict(g)
    if g["gate"] == "u3":
        if g["theta"] != "0*pi" or g["phi"] != "0*pi":
            raise ValueError("inverse rule only covers phase-only U3")
        g["lam"] = "pi/4" if g["lam"] == "-pi/4" else "-pi/4"
    elif g["gate"] not in {"h", "cx"}:
        raise ValueError("unexpected gate in RCCX inverse")
    return g


def rccx(a, b, t, adjoint):
    gates = [one("h", t), phase(t, 1), cx(b, t), phase(t, -1),
             cx(a, t), phase(t, 1), cx(b, t), phase(t, -1), one("h", t)]
    return [inverse(g) for g in reversed(gates)] if adjoint else gates


def variants():
    words = []
    for ps, pa, cs, ca in itertools.product((False, True), repeat=4):
        gates = [cx(0, 2), one("h", 0), cx(1, 3), one("h", 1)]
        gates += [cx(0, 1), cx(2, 3)]
        gates += [one("h", 2), cx(0, 2), one("h", 2)]
        a, b = (0, 3) if ps else (3, 0)
        gates += rccx(a, b, 2, pa)
        gates += [one("ry", 0, theta="-pi/4"), one("h", 0), cx(1, 0),
                  one("h", 0), one("ry", 0, theta="pi/4")]
        gates += [one("x", 1), one("ry", 2, theta="pi/4"), one("s", 2)]
        a, b = (3, 1) if cs else (1, 3)
        gates += rccx(a, b, 2, ca)
        gates += [one("sdg", 2), one("ry", 2, theta="-pi/4"), one("x", 1)]
        assert sum(g["gate"] == "cx" for g in gates) == 12
        vid = f"p_swap{int(ps)}_p_adj{int(pa)}_cch_swap{int(cs)}_cch_adj{int(ca)}"
        words.append({"n": 4, "gates": gates, "metadata": {
            "variant_id": vid, "choices": {"p_swap": ps, "p_adjoint": pa,
                                            "cch_swap": cs, "cch_adjoint": ca},
            "evidence_scope": SCOPE}})
    return words


def versions():
    return {"python": ".".join(map(str, sys.version_info[:3])), "numpy": np.__version__}


def dependencies():
    return {"numpy_init": digest(np.__file__),
            "numpy_multiarray_umath": digest(np._core._multiarray_umath.__file__)}


def source_hashes(root):
    return {"check_circuit.py": digest(root / "check_circuit.py")}


def prepare():
    if HERE.name != "numerical_sol" or HERE.parent.name != "work":
        raise RuntimeError("preparation only in collaboration/work/numerical_sol")
    root = HERE.parents[2]
    if any((HERE / name).exists() for name in ("config.json", "variants.json")):
        raise RuntimeError("prepared artifacts already exist; no implicit overwrite")
    dump(HERE / "variants.json", variants())
    config = dict(FIXED, script_sha256=digest(__file__), source_sha256=source_hashes(root),
                  runtime_versions=versions(), dependency_sha256=dependencies(),
                  variants_sha256=digest(HERE / "variants.json"))
    dump(HERE / "config.json", config)
    print(json.dumps({"prepared": True, "config_sha256": digest(HERE / "config.json"),
                      "script_sha256": config["script_sha256"], "variants_sha256": config["variants_sha256"],
                      "matrix_calls": 0}, sort_keys=True))


def guard():
    if not HERE.name.isdigit() or HERE.parent.name != "runs" or HERE.parent.parent.name != "experiments":
        raise RuntimeError("execute only from experiments/runs/<reserved-id>")
    root, run_id = HERE.parents[2], int(HERE.name)
    config = json.loads((HERE / "config.json").read_text())
    if set(config) != REQUIRED:
        raise RuntimeError(f"config keys differ: {sorted(set(config) ^ REQUIRED)}")
    for key, expected in FIXED.items():
        if config[key] != expected or type(config[key]) is not type(expected):
            raise RuntimeError(f"frozen config value mismatch: {key}")
    if config["script_sha256"] != digest(__file__):
        raise RuntimeError("script hash mismatch")
    if config["source_sha256"] != source_hashes(root):
        raise RuntimeError("checker source hash mismatch")
    if config["runtime_versions"] != versions() or config["dependency_sha256"] != dependencies():
        raise RuntimeError("runtime version or dependency hash mismatch")
    if config["variants_sha256"] != digest(HERE / "variants.json"):
        raise RuntimeError("variant manifest hash mismatch")
    words = json.loads((HERE / "variants.json").read_text())
    if words != variants() or len(words) != config["variant_count"]:
        raise RuntimeError("literal words do not match frozen generator")
    with sqlite3.connect(f"file:{root / 'experiments.sqlite3'}?mode=ro", uri=True) as db:
        row = db.execute("SELECT status,code_sha256,config_json,method FROM attempts WHERE id=?", (run_id,)).fetchone()
    if row is None or row[0] != "running" or row[1] != config["script_sha256"] or json.loads(row[2]) != config or row[3] != config["method"]:
        raise RuntimeError("ledger reservation does not match frozen script/config/method")
    with sqlite3.connect(f"file:{root / 'collaboration/board.sqlite3'}?mode=ro", uri=True) as db:
        board = db.execute("SELECT status,owner,run_id FROM hypotheses WHERE id=?", (config["board_id"],)).fetchone()
    if board != ("running", "numerical_sol", run_id):
        raise RuntimeError("active board ownership/run link mismatch")
    if any(os.environ[key] != "1" for key in THREAD_KEYS):
        raise RuntimeError("thread environment mismatch")
    return root, run_id, config, words


def evaluate(candidate, checker, tolerance, representatives):
    checked = checker.evaluate(candidate, tolerance=tolerance)
    U = np.eye(16, dtype=complex)
    for gate in candidate["gates"]:
        M = checker.cnot(gate["control"], gate["target"], 4) if gate["gate"] == "cx" else checker.embedded_one_qubit(checker.one_qubit_matrix(gate), gate["qubit"], 4)
        U = M @ U
    A = U @ checker.right_shift(4) @ U.conj().T
    diagonal = np.diag(A)
    off = A - np.diag(diagonal)
    roots = np.array([1, -1, 1j, -1j])
    distances = np.abs(diagonal[:, None] - roots[None, :])
    labels = np.argmin(distances, axis=1)
    counts = [int(np.count_nonzero(labels == k)) for k in range(4)]
    vid = candidate["metadata"]["variant_id"]
    duplicate = None
    for rid, R in representatives:
        if np.max(np.abs(U - R)) <= tolerance:
            duplicate = rid
            break
    if duplicate is None:
        representatives.append((vid, U.copy()))
    valid = checked["valid_diagonalizer"] and counts == FIXED["expected_multiplicities"]
    # UV=DU is equivalent for unitary U, but recorded directly as extra evidence.
    direct = U @ checker.right_shift(4) - np.diag(diagonal) @ U
    return {"variant_id": vid, "cx_count": checked["cnot_count"],
            "loss": float(np.vdot(off, off).real / 16), "maxoff": float(np.max(np.abs(off))),
            "uv_du_max_residual": float(np.max(np.abs(direct))),
            "diagonal_root_distance": float(np.max(np.min(distances, axis=1))),
            "nearest_diagonal_root_multiplicities": counts,
            "multiplicities_meaningful": bool(checked["valid_diagonalizer"]),
            "valid_diagonalizer": bool(valid), "same_unitary_as": duplicate,
            "shared_checker": checked}


def execute():
    started = time.monotonic()
    root, run_id, config, words = guard()
    # Exclusive marker is retained even on partial failure: no silent retries.
    with (HERE / "execution_started.json").open("x") as handle:
        json.dump({"run_id": run_id, "config_sha256": digest(HERE / "config.json"),
                   "utc_start": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, handle)
    remaining = config["time_limit_seconds"] - (time.monotonic() - started)
    if remaining <= 0:
        raise TimeoutError("preflight exhausted total budget")
    def timeout_handler(signum, frame):
        raise TimeoutError("frozen 30-second screen budget exceeded")
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, remaining)
    records, representatives = [], []
    try:
        sys.path.insert(0, str(root))
        import check_circuit as checker
        if Path(checker.__file__).resolve() != (root / "check_circuit.py").resolve():
            raise RuntimeError("unexpected loaded checker module path")
        out = HERE / "variants"
        out.mkdir(exist_ok=False)
        paths = []
        # Save all complete frozen words before the first matrix evaluation.
        for index, candidate in enumerate(words):
            path = out / f"{index:02d}_{candidate['metadata']['variant_id']}.json"
            dump(path, candidate)
            paths.append(path)
        for candidate, path in zip(words, paths):
            record = evaluate(candidate, checker, config["numeric_tolerance"], representatives)
            record.update(candidate_file=str(path.relative_to(HERE)), candidate_sha256=digest(path))
            records.append(record)
            dump(HERE / "partial_result.json", {"run_id": run_id, "records": records})
        best_index = min(range(len(records)), key=lambda i: (records[i]["loss"], records[i]["variant_id"]))
        dump(HERE / "best_candidate.json", words[best_index])
        dump(HERE / "all_metrics.json", records)
        result = {"status": "completed", "run_id": run_id, "config_sha256": digest(HERE / "config.json"),
                  "script_sha256": config["script_sha256"], "source_sha256": config["source_sha256"],
                  "dependency_sha256": config["dependency_sha256"], "runtime_versions": versions(),
                  "runtime_module_paths": {"numpy": np.__file__, "numpy_multiarray_umath": np._core._multiarray_umath.__file__, "check_circuit": checker.__file__},
                  "thread_environment": {key: os.environ[key] for key in THREAD_KEYS},
                  "evaluated_count": len(records), "distinct_unitaries_tolerance": config["numeric_tolerance"],
                  "distinct_unitary_count": len(representatives), "pass_count": sum(r["valid_diagonalizer"] for r in records),
                  "seconds": time.monotonic() - started, "best": records[best_index],
                  "best_candidate_sha256": digest(HERE / "best_candidate.json"),
                  "records": records, "evidence_scope": config["evidence_scope"]}
        dump(HERE / "result.json", result)
        print(json.dumps({key: result[key] for key in ("status", "run_id", "evaluated_count", "distinct_unitary_count", "pass_count", "seconds", "best")}, sort_keys=True))
    except BaseException as error:
        dump(HERE / "failure.json", {"run_id": run_id, "error_type": type(error).__name__,
                                     "error": str(error), "evaluated_count": len(records),
                                     "seconds": time.monotonic() - started})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true", help="serialize literals and freeze hashes; no matrices")
    args = parser.parse_args()
    prepare() if args.prepare else execute()
