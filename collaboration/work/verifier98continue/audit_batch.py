"""Independent finite preflight/postflight for the reserved 32-start full98 batch.

Uses the earlier verifier's NumPy/SciPy dense gate evaluator without changing it.
This script never runs an optimizer, computes derivatives, or perturbs a point.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE_PATH = ROOT / "collaboration/work/verifier98fresh/independent_audit.py"
BASE_SPEC = importlib.util.spec_from_file_location("verifier98fresh_base", BASE_PATH)
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)

EXPECTED_SEEDS = list(range(9261600, 9261632))
EXPECTED_SIGMAS = [0.3, 1.0, 2.0, 3.0]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_family(cfg: dict):
    spec = importlib.util.spec_from_file_location("full98_batch_audit_family", ROOT / cfg["family_module"])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def start_vector(seed: int, sigma: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.concatenate((rng.normal(0.0, sigma, 90), rng.uniform(-np.pi/2, np.pi/2, 8)))


def validate_plan(folder: Path, cfg: dict) -> dict:
    script = folder / "batch.py"
    code_hash = digest(script)
    if code_hash != cfg["code_sha256"]:
        raise SystemExit("batch source differs from frozen code hash")
    for relative, expected in cfg["dependencies"].items():
        if digest(ROOT / relative) != expected:
            raise SystemExit(f"frozen dependency changed: {relative}")
    if cfg["seeds"] != EXPECTED_SEEDS:
        raise SystemExit("seed list differs from 9261600..9261631")
    if cfg["local_sigmas"] != EXPECTED_SIGMAS:
        raise SystemExit("local sigma cycle differs from [0.3, 1, 2, 3]")
    limits = {
        "maxiter": 800, "maxfun": 20000,
        "internal_seconds": 285, "external_seconds": 300,
        "threads": 1, "nearhit_loss": 1e-16,
    }
    for key, expected in limits.items():
        if cfg.get(key) != expected:
            raise SystemExit(f"batch setting {key} differs from the frozen assignment")
    if cfg.get("family_module") != "collaboration/work/numerical98/hessian98.py":
        raise SystemExit("batch does not use the reviewed full98 coordinate family")
    source = script.read_text()
    for frozen_fragment in (
        "rng.normal(0,sigma,90)", "rng.uniform(-np.pi/2,np.pi/2,8)",
        "'maxls':30", "'ftol':1e-16", "'gtol':1e-12",
        "method='L-BFGS-B'",
    ):
        if frozen_fragment not in source:
            raise SystemExit(f"batch source does not implement frozen method detail: {frozen_fragment}")
    if len(cfg["seeds"]) != 32 or cfg.get("parameter_order", "").find("all coordinates free") < 0:
        raise SystemExit("batch must contain 32 starts with all 98 coordinates free")
    if (folder / "result.json").exists():
        raise SystemExit("batch already ran; refusing to label a post-run audit as preflight")
    return {
        "batch_script_sha256": code_hash,
        "config_sha256": digest(folder / "config.json"),
        "seeds": cfg["seeds"],
        "sigma_by_seed": {str(seed): EXPECTED_SIGMAS[i % 4] for i, seed in enumerate(EXPECTED_SEEDS)},
        "sampling_rule": "default_rng(seed): normal(0,sigma,90), then uniform(-pi/2,pi/2,8)",
        "limits": limits,
        "family_module": cfg["family_module"],
    }


def preflight(folder: Path, cfg: dict) -> dict:
    plan = validate_plan(folder, cfg)
    family = import_family(cfg)
    starts = []
    all_initial_points = []
    # Check reproducibility of all 32 exact points, while limiting independent
    # dense serialization/matrix checks to one representative for each sigma.
    for index, seed in enumerate(EXPECTED_SEEDS):
        sigma = EXPECTED_SIGMAS[index % 4]
        point = start_vector(seed, sigma)
        if point.shape != (98,) or not np.isfinite(point).all():
            raise SystemExit(f"invalid initial point for seed {seed}")
        if np.linalg.norm(point[84:87]) == 0 or np.linalg.norm(point[87:90]) == 0:
            raise SystemExit(f"seed {seed} has identity A1 or B1 in its sampled start")
        all_initial_points.append({
            "seed": seed, "sigma": sigma,
            "initial_vector_sha256_f64le": hashlib.sha256(np.asarray(point, dtype="<f8").tobytes()).hexdigest(),
            "A1_norm": float(np.linalg.norm(point[84:87])),
            "B1_norm": float(np.linalg.norm(point[87:90])),
        })
        if index not in range(4):
            continue
        native = family.serialize(point)
        compiled = family.helper.compile_native(native)
        coordinate = base.coordinate_matrix(point)
        nmat, ng = base.gate_list(native)
        cmat, cg = base.gate_list(compiled)
        matrix_err = base.phase_error(coordinate, nmat)
        compile_err = base.phase_error(nmat, cmat)
        if matrix_err > 2e-12 or compile_err > 2e-12:
            raise SystemExit(f"prefit circuit mismatch for seed {seed}: {matrix_err}, {compile_err}")
        if (sum(g["gate"] == "cx" for g in ng),
                sum(g["gate"] == "xx_yy" for g in ng),
                sum(g["gate"] == "cx" for g in cg)) != (4, 4, 12):
            raise SystemExit(f"gate count mismatch for seed {seed}")
        starts.append({
            "seed": seed, "sigma": sigma,
            "initial_vector_sha256_f64le": hashlib.sha256(np.asarray(point, dtype="<f8").tobytes()).hexdigest(),
            "A1_norm": float(np.linalg.norm(point[84:87])),
            "B1_norm": float(np.linalg.norm(point[87:90])),
            "native_cx_count": 4, "native_xx_yy_count": 4, "compiled_cx_count": 12,
            "coordinate_vs_native_matrix_error_up_to_phase": matrix_err,
            "native_vs_compiled_matrix_error_up_to_phase": compile_err,
            "native_unitarity_max_error": float(np.max(np.abs(nmat.conj().T @ nmat - np.eye(16)))),
        })
    plan["prefit_starts"] = starts
    plan["all_32_initial_points_reproducible"] = True
    plan["all_32_initial_point_hashes_and_interior_norms"] = all_initial_points
    plan["dense_matrix_checks"] = "one seed per each of four sigma scales"
    plan["run_id"] = int(folder.name)
    plan["board_hypothesis_id"] = 78
    plan["parameter_order"] = cfg["parameter_order"]
    plan["native_schedule"] = cfg["native_schedule"]
    plan["prefit_status"] = "APPROVED for exactly one execution of this reserved script/config; threshold is a promising trigger, not exact success"
    plan["scope"] = "32 deterministic initial points and serialization checks only; no objective, fit, derivative, or perturbation"
    return plan


def postfit(folder: Path, cfg: dict, approval: dict) -> dict:
    result_path = folder / "result.json"
    result = json.loads(result_path.read_text())
    if result.get("run_id") != int(folder.name):
        raise SystemExit("batch result run ID differs from reserved workspace")
    if result.get("config") != cfg:
        raise SystemExit("result configuration differs from frozen workspace config")
    rows = result.get("rows")
    if not isinstance(rows, list) or not rows:
        raise SystemExit("batch result contains no completed seed endpoints")
    seeds = [row.get("seed") for row in rows]
    if seeds != EXPECTED_SEEDS[:len(seeds)]:
        raise SystemExit("completed seed rows are not a prefix of the frozen sequence")
    if len(rows) > len(EXPECTED_SEEDS):
        raise SystemExit("batch contains more than 32 seed endpoints")
    batch_elapsed = float(result.get("elapsed_seconds", float("inf")))
    if batch_elapsed > cfg["internal_seconds"] or batch_elapsed > cfg["external_seconds"]:
        raise SystemExit("recorded batch elapsed time exceeds a frozen wall limit")
    preflight_hashes = {item["seed"]: item["initial_vector_sha256_f64le"]
                        for item in approval["all_32_initial_point_hashes_and_interior_norms"]}
    expected_names = set()
    audited = []
    for row_index, row in enumerate(rows):
        seed = EXPECTED_SEEDS[row_index]
        sigma = EXPECTED_SIGMAS[row_index % 4]
        x0 = start_vector(seed, sigma)
        actual_initial = np.asarray(row["initial"], dtype=float)
        if actual_initial.shape != (98,) or not np.array_equal(actual_initial, x0):
            raise SystemExit(f"seed {seed} initial vector does not reproduce from frozen RNG rule")
        actual_initial_hash = hashlib.sha256(np.asarray(actual_initial, dtype="<f8").tobytes()).hexdigest()
        if actual_initial_hash != preflight_hashes.get(seed):
            raise SystemExit(f"seed {seed} initial vector differs from the preflight-frozen point")
        if row.get("sigma") != sigma:
            raise SystemExit(f"seed {seed} sigma field differs from cyclic schedule")
        if int(row["evaluations"]) > cfg["maxfun"]:
            raise SystemExit(f"seed {seed} exceeded maxfun")
        if row.get("nit") is not None and int(row["nit"]) > cfg["maxiter"]:
            raise SystemExit(f"seed {seed} exceeded maxiter")
        endpoint = np.asarray(row["best"], dtype=float)
        if endpoint.shape != (98,) or not np.isfinite(endpoint).all():
            raise SystemExit(f"seed {seed} saved endpoint is not a finite 98-vector")
        native_path = folder / row["native"]
        compiled_path = folder / row["compiled"]
        expected_names.update((native_path.name, compiled_path.name))
        if digest(native_path) != row["native_sha256"] or digest(compiled_path) != row["compiled_sha256"]:
            raise SystemExit(f"seed {seed} saved gate-list hash mismatch")
        nmat, ng = base.gate_list(native_path)
        cmat, cg = base.gate_list(compiled_path)
        coordinate = base.coordinate_matrix(endpoint)
        coordinate_err = base.phase_error(coordinate, nmat)
        compile_err = base.phase_error(nmat, cmat)
        unitarity_n = float(np.max(np.abs(nmat.conj().T @ nmat - np.eye(16))))
        unitarity_c = float(np.max(np.abs(cmat.conj().T @ cmat - np.eye(16))))
        target = nmat @ base.right_shift() @ nmat.conj().T
        off = target - np.diag(np.diag(target))
        fro = float(np.linalg.norm(off))
        maxoff = float(np.max(np.abs(off)))
        loss = fro * fro / 16
        counts = (sum(g["gate"] == "cx" for g in ng),
                  sum(g["gate"] == "xx_yy" for g in ng),
                  sum(g["gate"] == "cx" for g in cg))
        if counts != (4, 4, 12):
            raise SystemExit(f"seed {seed} circuit counts {counts}, expected (4,4,12)")
        if max(coordinate_err, compile_err, unitarity_n, unitarity_c) > 2e-12:
            raise SystemExit(f"seed {seed} matrix validation failed")
        if abs(loss - float(row["loss"])) > 2e-12:
            raise SystemExit(f"seed {seed} independently recomputed loss disagrees with saved objective")
        audited.append({
            "seed": seed, "sigma": sigma, "loss_saved": float(row["loss"]),
            "loss_independent": loss, "loss_abs_difference": abs(loss-float(row["loss"])),
            "gradient_norm_saved": float(row["gradient_norm"]),
            "evaluations": int(row["evaluations"]), "iterations": row.get("nit"),
            "termination": row["termination"], "target_offdiagonal_frobenius": fro,
            "target_offdiagonal_max_entry": maxoff,
            "coordinate_vs_native_matrix_error_up_to_phase": coordinate_err,
            "native_vs_compiled_matrix_error_up_to_phase": compile_err,
            "native_unitarity_max_error": unitarity_n,
            "compiled_unitarity_max_error": unitarity_c,
            "native_cx_count": counts[0], "native_xx_yy_count": counts[1],
            "compiled_cx_count": counts[2],
            "native_sha256": digest(native_path), "compiled_sha256": digest(compiled_path),
        })
    actual_endpoint_files = {p.name for p in folder.glob("endpoint_*_native.json")} | {p.name for p in folder.glob("endpoint_*_compiled.json")}
    if actual_endpoint_files != expected_names:
        raise SystemExit(f"unmatched endpoint files: expected {sorted(expected_names)}, found {sorted(actual_endpoint_files)}")
    best_index = min(range(len(audited)), key=lambda i: audited[i]["loss_independent"])
    if result.get("best_index") != best_index:
        raise SystemExit("saved best_index differs from independent objective ranking")
    best_row = rows[best_index]
    for kind in ("native", "compiled"):
        global_path = folder / f"best_{kind}.json"
        if not global_path.exists() or global_path.read_bytes() != (folder / best_row[kind]).read_bytes():
            raise SystemExit(f"global best_{kind} list is not the minimum-loss row's endpoint")
    early_nearhit = result.get("stop_reason") == "promising_nearhit_requires_independent_verification"
    if early_nearhit and audited[-1]["loss_independent"] >= cfg["nearhit_loss"]:
        raise SystemExit("batch claims nearhit stop without endpoint meeting threshold")
    if not early_nearhit and any(x["loss_independent"] < cfg["nearhit_loss"] for x in audited):
        raise SystemExit("batch failed to stop at its frozen promising-nearhit threshold")
    approval.update({
        "run_id": int(folder.name), "result_sha256": digest(result_path),
        "completed_starts": len(rows), "stop_reason": result.get("stop_reason"),
        "elapsed_seconds": batch_elapsed,
        "best_seed": audited[best_index]["seed"], "best_loss": audited[best_index]["loss_independent"],
        "any_numerically_passing_endpoint": any(x["target_offdiagonal_max_entry"] < 1e-10 for x in audited),
        "nearhit_threshold": cfg["nearhit_loss"], "starts": audited,
        "scope": "independent matrix/objective verification of every saved per-seed best endpoint; numerical checks are not exact certification",
    })
    return approval


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: python3 audit_batch.py RESERVED_RUN_DIR [postfit]")
    folder = Path(sys.argv[1]).resolve()
    cfg = json.loads((folder / "config.json").read_text())
    if len(sys.argv) == 2:
        report = preflight(folder, cfg)
        out = Path(__file__).resolve().parent / "PREFIT_APPROVAL.json"
    elif sys.argv[2] == "postfit":
        approval_path = Path(__file__).resolve().parent / "PREFIT_APPROVAL.json"
        report = postfit(folder, cfg, json.loads(approval_path.read_text()))
        out = Path(__file__).resolve().parent / "POSTFIT_AUDIT.json"
    else:
        raise SystemExit("phase must be 'postfit'")
    report["auditor_script_sha256"] = digest(Path(__file__))
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in report if k != "prefit_starts" and k != "starts"}, indent=2))


if __name__ == "__main__":
    main()
