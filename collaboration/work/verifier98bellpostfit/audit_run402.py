"""Independent finite audit of run402 endpoints; no fit or derivative calls."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "experiments/runs/402"
PREFLIGHT_PATH = ROOT / "experiments/runs/404/preflight.py"
EXPECTED_PREFLIGHT_SHA = "9422a752053361ede634b2dde88933def820c7525957724412684e33975cca92"
EXPECTED_CONFIG_SHA = "1fb2864e4f9750324080d95f8193d46dce75c3939373036b4187b8e080c4616a"
EXPECTED_FAMILY_SHA = "9936da7408d5202a00e784daa4b33d950efa68229bec06c006457862411303c8"
EXPECTED_SEARCH_SHA = "3f4777c2492e34e47fb5655df4c288f432f8dbad0e74230afeff872c257448ac"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def import_independent():
    if sha(PREFLIGHT_PATH) != EXPECTED_PREFLIGHT_SHA:
        raise SystemExit("independent preflight helper source hash changed")
    spec = importlib.util.spec_from_file_location("audited_independent_numpy", PREFLIGHT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def v4():
    out = np.zeros((16, 16), dtype=complex)
    for col in range(16):
        out[((col & 1) << 3) | (col >> 1), col] = 1.0
    return out


def loss_metrics(unitary, target):
    transformed = unitary @ target @ unitary.conj().T
    diagonal = np.diag(np.diag(transformed))
    off = transformed - diagonal
    return {
        "normalized_offdiagonal_frobenius_loss": float(np.linalg.norm(off, "fro") ** 2 / 16),
        "offdiagonal_frobenius_norm": float(np.linalg.norm(off, "fro")),
        "offdiagonal_max_entry": float(np.max(np.abs(off))),
        "unitarity_max_error": float(np.max(np.abs(unitary.conj().T @ unitary - np.eye(16)))),
    }


def validate_native_word(gates, native_word):
    cursor = 0
    def layer(label):
        nonlocal cursor
        for q in range(4):
            if cursor >= len(gates) or gates[cursor].get("gate") != "u3" or gates[cursor].get("qubit") != q:
                raise ValueError(f"{label} not serialized in chronological q0..q3 order")
            cursor += 1
    def entangler(kind, a, b):
        nonlocal cursor
        if cursor >= len(gates): raise ValueError("native list ended before declared word")
        g = gates[cursor]
        ok = ((kind == "cx" and g.get("gate") == "cx" and g.get("control") == a and g.get("target") == b)
              or (kind == "xx_yy" and g.get("gate") == "xx_yy" and g.get("qubits") == [a,b]))
        if not ok: raise ValueError(f"native entangler mismatch at position {cursor}")
        cursor += 1
    def local_special(label, q):
        nonlocal cursor
        if cursor >= len(gates) or gates[cursor].get("gate") != "u3" or gates[cursor].get("qubit") != q:
            raise ValueError(f"{label} not serialized at declared chronology")
        cursor += 1
    for tok in native_word:
        if tok[0] == "layer": layer(tok[1])
        elif tok[0] in ("cx", "xx_yy"): entangler(tok[0],tok[1],tok[2])
        elif tok[0] in ("A1", "B1"): local_special(tok[0],tok[1])
        else: raise ValueError(f"unknown word token: {tok}")
    if cursor != len(gates): raise ValueError(f"{len(gates)-cursor} undeclared native gate(s)")


def audit_point(name, vector, native_path, compiled_path, result, cfg, independent):
    vector = np.asarray(vector, dtype=np.float64)
    if vector.shape != (98,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} vector is not a finite 98-vector")
    native_sha, compiled_sha = sha(native_path), sha(compiled_path)
    native_data = json.loads(Path(native_path).read_text())
    compiled_data = json.loads(Path(compiled_path).read_text())
    if native_sha != result[f"{name}_native_sha256"] or compiled_sha != result[f"{name}_compiled_sha256"]:
        raise ValueError(f"{name} serialized list hash differs from run result")
    if native_data.get("n") != 4 or compiled_data.get("n") != 4:
        raise ValueError(f"{name} gate list does not declare four qubits")
    validate_native_word(native_data["gates"], cfg["native_word"])
    native_cx = sum(g.get("gate") == "cx" for g in native_data["gates"])
    native_f = sum(g.get("gate") == "xx_yy" for g in native_data["gates"])
    compiled_cx = sum(g.get("gate") == "cx" for g in compiled_data["gates"])
    if native_cx != 4 or native_f != 4 or compiled_cx != 12:
        raise ValueError(f"{name} native/compiled entangler count mismatch")
    if any(g.get("gate") == "xx_yy" for g in compiled_data["gates"]):
        raise ValueError(f"{name} compiled list still contains native XX/YY")

    coordinate = independent.coordinate_matrix(vector)
    native = independent.gate_list_matrix(native_data)
    compiled = independent.gate_list_matrix(compiled_data)
    target = v4()
    metrics = {
        "coordinate": loss_metrics(coordinate, target),
        "native": loss_metrics(native, target),
        "compiled": loss_metrics(compiled, target),
        "coordinate_vs_native_phase_error": independent.phase_error(coordinate, native),
        "coordinate_vs_compiled_phase_error": independent.phase_error(coordinate, compiled),
        "native_vs_compiled_phase_error": independent.phase_error(native, compiled),
        "native_sha256": native_sha,
        "compiled_sha256": compiled_sha,
        "native_gate_count": len(native_data["gates"]),
        "compiled_gate_count": len(compiled_data["gates"]),
        "native_cx_count": native_cx,
        "native_xx_yy_count": native_f,
        "compiled_cx_count": compiled_cx,
    }
    # The run's in-process checker is supporting provenance, not the independent result.
    checker_key = f"{name}_native_checker"
    checker = result.get(checker_key)
    if checker is not None:
        metrics["run_checker_loss"] = checker.get("normalized_loss")
        metrics["run_checker_offdiagonal_max"] = checker.get("off_diagonal_error")
    if max(metrics[k] for k in ("coordinate_vs_native_phase_error", "coordinate_vs_compiled_phase_error", "native_vs_compiled_phase_error")) > 2e-10:
        raise ValueError(f"{name} saved list does not match its complete coordinate vector")
    return metrics


def main():
    if sha(RUN / "config.json") != EXPECTED_CONFIG_SHA or sha(RUN / "family.py") != EXPECTED_FAMILY_SHA or sha(RUN / "search.py") != EXPECTED_SEARCH_SHA:
        raise SystemExit("run402 frozen fit inputs changed")
    cfg = json.loads((RUN / "config.json").read_text())
    result = json.loads((RUN / "result.json").read_text())
    if result.get("run_id") != 402 or result.get("config") != cfg:
        raise SystemExit("run result/config identity mismatch")
    independent = import_independent()
    initial = np.asarray(result["initial98"], dtype=np.float64)
    base = independent.expected_base()
    noise = np.random.default_rng(cfg["seed"]).normal(0.0, cfg["noise_sigma"], size=98)
    expected_initial = base + noise
    initial_hash = hashlib.sha256(np.asarray(initial, dtype="<f8").tobytes()).hexdigest()
    preflight = json.loads((ROOT / "experiments/runs/404/preflight_result.json").read_text())
    preflight_initial = next(p for p in preflight["points"] if p["name"] == "seeded_initial")
    if not np.array_equal(initial, expected_initial) or initial_hash != preflight_initial["point_sha256_f64le"]:
        raise SystemExit("run402 initial vector differs from frozen RNG start/preflight point")
    if result.get("seed") != cfg["seed"] or len(result.get("gaussian_perturbation98", [])) != 98:
        raise SystemExit("seeded initialization record malformed")
    if not np.array_equal(np.asarray(result["gaussian_perturbation98"], dtype=float), noise):
        raise SystemExit("saved perturbation differs from seeded NumPy vector")
    if not result.get("best_evaluated") or len(result.get("best98", [])) != 98:
        raise SystemExit("run has no complete evaluated best point")

    initial_data = audit_point("initial", initial, RUN/"initial_native.json", RUN/"initial_compiled.json", result, cfg, independent)
    best_data = audit_point("best", np.asarray(result["best98"], dtype=float), RUN/"best_native.json", RUN/"best_compiled.json", result, cfg, independent)
    evaluations = result["evaluations"]
    history = result["history"]
    if evaluations != len(history) or evaluations > cfg["maxfun"]:
        raise SystemExit("optimizer evaluation history/count exceeds frozen budget")
    if [h["evaluation"] for h in history] != list(range(1, evaluations + 1)):
        raise SystemExit("optimizer history evaluation indices are not contiguous")
    min_record = min(history, key=lambda h: h["loss"])
    if result.get("best_evaluation") != min_record["evaluation"] or abs(result["best_loss"] - min_record["loss"]) > 1e-12:
        raise SystemExit("recorded best point does not match minimum history loss")
    if result.get("optimizer", {}).get("nfev") != evaluations:
        raise SystemExit("SciPy nfev differs from guarded objective call count")
    if result["optimizer"].get("nit", 10**9) > cfg["maxiter"] or result.get("elapsed_seconds", 10**9) > cfg["internal_seconds"]:
        raise SystemExit("fit exceeded frozen iteration/internal-time budget")
    if result.get("final_recomputation_count_outside_optimizer_maxfun") not in (0, 1):
        raise SystemExit("outside-maxfun final recomputation count exceeds one")
    if abs(best_data["coordinate"]["normalized_offdiagonal_frobenius_loss"] - result["best_loss_recomputed"]) > 2e-10:
        raise SystemExit("independent best loss disagrees with saved final recomputation")
    output = {
        "run_id": 402,
        "board_hypothesis_id": 108,
        "schema": "independent-reordered-full98-postfit-audit-v1",
        "fit_source_hashes": {"config": EXPECTED_CONFIG_SHA, "family": EXPECTED_FAMILY_SHA, "search": EXPECTED_SEARCH_SHA},
        "preflight_result_sha256": sha(ROOT / "experiments/runs/404/preflight_result.json"),
        "initial_vector_sha256_f64le": initial_hash,
        "evaluations": evaluations,
        "history_count": len(history),
        "best_evaluation": result["best_evaluation"],
        "best_loss_recorded": result["best_loss"],
        "best_loss_recomputed": result["best_loss_recomputed"],
        "best_gradient_norm_recorded": result["best_gradient_norm"],
        "fit_elapsed_seconds": result["elapsed_seconds"],
        "fit_termination": result["termination"],
        "initial_endpoint": initial_data,
        "best_endpoint": best_data,
        "budget": {"maxiter": cfg["maxiter"], "maxfun": cfg["maxfun"], "iterations": result["optimizer"]["nit"], "evaluations": evaluations,
                   "internal_seconds": cfg["internal_seconds"], "elapsed_seconds": result["elapsed_seconds"], "threads": cfg["threads"]},
        "scope": "Independent finite validation of initial/best saved full gate lists and histories only; no optimizer, AD, derivative, assignment or fit call.",
    }
    own_out = Path(__file__).resolve().parent / "POSTFIT402_AUDIT.json"
    run_out = RUN / "POSTFIT402_AUDIT.json"
    if own_out.exists() or run_out.exists():
        raise SystemExit("refusing to repeat post-fit endpoint audit")
    for out in (own_out, run_out):
        out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
