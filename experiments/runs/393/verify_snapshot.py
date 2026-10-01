"""Independent NumPy/SciPy audit for the fresh full98 fit.

This module performs finite matrix evaluation only. It never differentiates,
optimizes, perturbs the start, or classifies a basin. Run it only after the
numerical researcher freezes a complete artifact directory.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[3]
PAIRS = ((0, 2), (1, 3), (0, 1), (2, 3))
CX_PAIRS = ((0, 1), (2, 1), (3, 1), (1, 2))
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = {"rx": X, "ry": Y, "rz": Z}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def phase_error(a: np.ndarray, b: np.ndarray) -> float:
    p = np.vdot(b, a)
    if abs(p) == 0:
        return float("inf")
    return float(np.max(np.abs(a - (p / abs(p)) * b)))


def su2(axis: np.ndarray) -> np.ndarray:
    """exp(-i axis.sigma/2), evaluated independently in NumPy."""
    x, y, z = map(float, axis)
    theta = math.sqrt(x*x + y*y + z*z)
    if theta == 0:
        return I2.copy()
    return (math.cos(theta/2) * I2
            - 1j * math.sin(theta/2) / theta * (x*X + y*Y + z*Z))


def embed_one(gate: np.ndarray, q: int) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for k in range(4):
        out = np.kron(out, gate if k == q else I2)
    return out


def cx(control: int, target: int) -> np.ndarray:
    out = np.zeros((16, 16), dtype=complex)
    for col in range(16):
        row = col
        if (col >> (3-control)) & 1:
            row ^= 1 << (3-target)
        out[row, col] = 1
    return out


def f_gate(a: float, b: float, q0: int, q1: int) -> np.ndarray:
    xx = embed_one(X, q0) @ embed_one(X, q1)
    yy = embed_one(Y, q0) @ embed_one(Y, q1)
    return expm(1j * (float(a) * xx + float(b) * yy))


def coordinate_matrix(point: np.ndarray) -> np.ndarray:
    """Evaluate the declared 98-coordinate chronological family directly."""
    z = np.asarray(point, dtype=float)
    if z.shape != (98,):
        raise ValueError(f"expected point shape (98,), got {z.shape}")
    local = z[:84].reshape(7, 4, 3)
    ab = z[84:90].reshape(2, 3)
    f = z[90:].reshape(4, 2)
    u = np.eye(16, dtype=complex)

    def layer(k: int) -> None:
        nonlocal u
        for q in range(4):
            u = embed_one(su2(local[k, q]), q) @ u

    layer(0)
    u = cx(*CX_PAIRS[0]) @ u
    u = embed_one(su2(ab[0]), 1) @ u
    u = cx(*CX_PAIRS[1]) @ u
    u = embed_one(su2(ab[1]), 1) @ u
    u = cx(*CX_PAIRS[2]) @ u
    layer(1)
    u = f_gate(*f[0], *PAIRS[0]) @ u
    layer(2)
    u = f_gate(*f[1], *PAIRS[1]) @ u
    layer(3)
    u = cx(*CX_PAIRS[3]) @ u
    layer(4)
    u = f_gate(*f[2], *PAIRS[2]) @ u
    layer(5)
    u = f_gate(*f[3], *PAIRS[3]) @ u
    layer(6)
    return u


def gate_matrix(g: dict) -> np.ndarray:
    kind = g["gate"]
    if kind == "cx":
        return cx(int(g["control"]), int(g["target"]))
    if kind == "xx_yy":
        return f_gate(float(g["a"]), float(g["b"]), *map(int, g["qubits"]))
    if kind == "u3":
        th, ph, la = map(float, (g["theta"], g["phi"], g["lam"]))
        one = np.array([
            [np.cos(th/2), -np.exp(1j*la)*np.sin(th/2)],
            [np.exp(1j*ph)*np.sin(th/2), np.exp(1j*(ph+la))*np.cos(th/2)],
        ], dtype=complex)
        return embed_one(one, int(g["qubit"]))
    if kind in PAULI:
        return embed_one(expm(-0.5j * float(g["theta"]) * PAULI[kind]), int(g["qubit"]))
    raise ValueError(f"unsupported gate {g}")


def gate_list(payload: dict | Path) -> tuple[np.ndarray, list[dict]]:
    data = payload if isinstance(payload, dict) else json.loads(payload.read_text())
    if data.get("n") != 4 or not isinstance(data.get("gates"), list):
        raise ValueError("gate list must declare n=4 and a gates array")
    u = np.eye(16, dtype=complex)
    for g in data["gates"]:
        u = gate_matrix(g) @ u
    return u, data["gates"]


def right_shift() -> np.ndarray:
    v = np.zeros((16, 16), dtype=complex)
    for col in range(16):
        v[((col & 1) << 3) | (col >> 1), col] = 1
    return v


def evaluate_lists(folder: Path, stem: str) -> dict:
    native_path = folder / f"{stem}_native.json"
    compiled_path = folder / f"{stem}_compiled.json"
    nmat, ng = gate_list(native_path)
    cmat, cg = gate_list(compiled_path)
    v = right_shift()
    b = nmat @ v @ nmat.conj().T
    off = b - np.diag(np.diag(b))
    return {
        "native_path": str(native_path.relative_to(ROOT)),
        "compiled_path": str(compiled_path.relative_to(ROOT)),
        "native_sha256": sha256(native_path),
        "compiled_sha256": sha256(compiled_path),
        "native_gate_count": len(ng),
        "native_cx_count": sum(g["gate"] == "cx" for g in ng),
        "native_xx_yy_count": sum(g["gate"] == "xx_yy" for g in ng),
        "compiled_gate_count": len(cg),
        "compiled_cx_count": sum(g["gate"] == "cx" for g in cg),
        "native_compiled_matrix_error_up_to_phase": phase_error(nmat, cmat),
        "native_unitarity_max_error": float(np.max(np.abs(nmat.conj().T @ nmat - np.eye(16)))),
        "compiled_unitarity_max_error": float(np.max(np.abs(cmat.conj().T @ cmat - np.eye(16)))),
        "target_offdiagonal_frobenius": float(np.linalg.norm(off)),
        "target_offdiagonal_max_entry": float(np.max(np.abs(off))),
        "normalized_loss": float(np.linalg.norm(off)**2 / 16),
    }


def static_review(folder: Path, cfg: dict) -> dict:
    script_path = folder / "search.py"
    if sha256(script_path) != cfg["code_sha256"]:
        raise SystemExit("numerical script hash does not match frozen config")
    for rel, digest in cfg["source_hashes"].items():
        source = folder / rel if rel == "config.json" else ROOT / rel
        if sha256(source) != digest:
            raise SystemExit(f"frozen dependency hash mismatch: {rel}")
    limits = {"maxiter": 300, "maxfun": 10000, "maxls": 30, "ftol": 1e-15,
              "gtol": 1e-10, "internal_timeout_seconds": 110,
              "external_timeout_seconds": 120, "threads": 1}
    for key, value in limits.items():
        if cfg.get(key) != value:
            raise SystemExit(f"fit limit {key} differs from the frozen hypothesis")
    if cfg.get("seed") != 9261501 or cfg.get("parameters") != 98:
        raise SystemExit("seed or parameter count differs from the frozen hypothesis")
    if cfg.get("starts") != 1 or cfg.get("restarts") != 0:
        raise SystemExit("fit is not frozen to one start with no restarts")
    point = np.asarray(cfg["initial_vector"], dtype=float)
    expected = np.random.default_rng(9261501).normal(0, 0.3, 98)
    if point.shape != (98,) or not np.array_equal(point, expected):
        raise SystemExit("start is not exact default_rng(9261501).normal(0,.3,98)")
    a_norm, b_norm = float(np.linalg.norm(point[84:87])), float(np.linalg.norm(point[87:90]))
    if a_norm <= 1e-12 or b_norm <= 1e-12:
        raise SystemExit("A1 or B1 is identity in the frozen start")
    return {
        "hypothesis": 76,
        "script_sha256": sha256(script_path),
        "config_sha256": sha256(folder / "config.json"),
        "seed": 9261501,
        "initial_vector_shape": list(point.shape),
        "initial_vector_sha256_f64le": hashlib.sha256(np.asarray(point, dtype="<f8").tobytes()).hexdigest(),
        "initial_rng_match": True,
        "A1_vector": point[84:87].tolist(), "A1_norm": a_norm,
        "B1_vector": point[87:90].tolist(), "B1_norm": b_norm,
        "all_98_coordinates_free": True,
        "budget": limits,
        "scope": "independent provenance review; no objective, fit, derivatives, or perturbations",
    }


def prefit(folder: Path, cfg: dict) -> dict:
    """Call after ledger reservation and before the optimizer is started."""
    module_path = ROOT / cfg["family_module"]
    spec = importlib.util.spec_from_file_location("fresh98_family_for_prefit", module_path)
    family = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(family)
    point = np.asarray(cfg["initial_vector"], dtype=float)
    native, compiled = family.serialize(point), family.helper.compile_native(family.serialize(point))
    coordinate_matrix_value = coordinate_matrix(point)
    native_matrix, native_gates = gate_list(native)
    compiled_matrix, compiled_gates = gate_list(compiled)
    report = static_review(folder, cfg)
    report.update({
        "native_gate_count": len(native_gates),
        "native_cx_count": sum(g["gate"] == "cx" for g in native_gates),
        "native_xx_yy_count": sum(g["gate"] == "xx_yy" for g in native_gates),
        "compiled_gate_count": len(compiled_gates),
        "compiled_cx_count": sum(g["gate"] == "cx" for g in compiled_gates),
        "coordinate_vs_native_matrix_error_up_to_phase": phase_error(coordinate_matrix_value, native_matrix),
        "native_vs_compiled_matrix_error_up_to_phase": phase_error(native_matrix, compiled_matrix),
        "native_unitarity_max_error": float(np.max(np.abs(native_matrix.conj().T @ native_matrix - np.eye(16)))),
        "native_schedule": cfg["native_schedule"],
        "chronology": cfg["chronology"],
        "scope": "prefit independent NumPy/SciPy initial-gate evaluation; no optimizer or derivatives",
    })
    if (report["native_cx_count"], report["native_xx_yy_count"], report["compiled_cx_count"]) != (4, 4, 12):
        raise SystemExit("initial gate counts differ from 4CX+4F and compiled12CX")
    if report["coordinate_vs_native_matrix_error_up_to_phase"] > 2e-12:
        raise SystemExit("coordinate matrix disagrees with serialized native initial circuit")
    if report["native_vs_compiled_matrix_error_up_to_phase"] > 2e-12:
        raise SystemExit("serialized native and compiled initial circuits disagree")
    return report


def postfit(folder: Path, cfg: dict, prior: dict) -> dict:
    result_path = folder / "result.json"
    result = json.loads(result_path.read_text())
    if result.get("run_id") != int(folder.name) or result.get("hypothesis_id") != 75:
        raise SystemExit("result does not match reserved run or hypothesis 75")
    if result.get("seed") != 9261501:
        raise SystemExit("result seed differs from prefit approval")
    counts = result.get("counts", {})
    if counts.get("starts") != 1 or counts.get("restarts") != 0 or counts.get("threads") != 1:
        raise SystemExit("result violates single-start/no-restart/single-thread budget")
    for label, vector_key in (("initial", "initial_vector"), ("best", "best_vector")):
        vector = np.asarray(result[vector_key], dtype=float)
        if vector.shape != (98,):
            raise SystemExit(f"{vector_key} must have 98 entries")
        if label == "initial" and not np.array_equal(vector, np.asarray(cfg["initial_vector"], dtype=float)):
            raise SystemExit("result initial vector differs from the frozen vector")
        npmat = coordinate_matrix(vector)
        npath, cpath = folder / f"{label}_native.json", folder / f"{label}_compiled.json"
        nmat, ng = gate_list(npath)
        cmat, cg = gate_list(cpath)
        if result.get(f"{label}_native") != json.loads(npath.read_text()):
            raise SystemExit(f"{label} native gate-list file differs from result.json payload")
        if result.get(f"{label}_compiled") != json.loads(cpath.read_text()):
            raise SystemExit(f"{label} compiled gate-list file differs from result.json payload")
        offmat = nmat @ right_shift() @ nmat.conj().T
        offmat -= np.diag(np.diag(offmat))
        audit = {
            "native_sha256": sha256(npath), "compiled_sha256": sha256(cpath),
            "native_gate_count": len(ng),
            "native_cx_count": sum(g["gate"] == "cx" for g in ng),
            "native_xx_yy_count": sum(g["gate"] == "xx_yy" for g in ng),
            "compiled_gate_count": len(cg),
            "compiled_cx_count": sum(g["gate"] == "cx" for g in cg),
            "coordinate_vs_native_matrix_error_up_to_phase": phase_error(npmat, nmat),
            "native_vs_compiled_matrix_error_up_to_phase": phase_error(nmat, cmat),
            "native_unitarity_max_error": float(np.max(np.abs(nmat.conj().T @ nmat - np.eye(16)))),
            "compiled_unitarity_max_error": float(np.max(np.abs(cmat.conj().T @ cmat - np.eye(16)))),
            "target_offdiagonal_frobenius": float(np.linalg.norm(offmat)),
            "target_offdiagonal_max_entry": float(np.max(np.abs(offmat))),
            "normalized_loss": float(np.linalg.norm(offmat)**2 / 16),
        }
        if (audit["native_cx_count"], audit["native_xx_yy_count"], audit["compiled_cx_count"]) != (4, 4, 12):
            raise SystemExit(f"{label} gate counts differ from 4CX+4F/compiled12CX")
        if audit["coordinate_vs_native_matrix_error_up_to_phase"] > 2e-12:
            raise SystemExit(f"{label} coordinates and native gate list disagree")
        if audit["native_vs_compiled_matrix_error_up_to_phase"] > 2e-12:
            raise SystemExit(f"{label} native and compiled lists disagree")
        prior[label] = audit
    prior["run_id"] = int(folder.name)
    prior["postfit_audit_script_sha256"] = sha256(Path(__file__))
    prior["fit_result_sha256"] = sha256(result_path)
    prior["saved_loss_initial_difference"] = abs(prior["initial"]["normalized_loss"] - result["initial_loss"])
    prior["saved_loss_best_difference"] = abs(prior["best"]["normalized_loss"] - result["best_loss"])
    prior["scope"] = "postfit independent complete-list validation; numerical evidence, not an exact certificate or family exclusion"
    return prior


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: python3 independent_audit.py ARTIFACT_DIR [postfit]")
    folder = Path(sys.argv[1]).resolve()
    cfg = json.loads((folder / "config.json").read_text())
    result = postfit(folder, cfg, json.loads((folder / "PREFIT_APPROVAL.json").read_text())) if len(sys.argv) == 3 else prefit(folder, cfg)
    out = folder / "independent_audit.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
