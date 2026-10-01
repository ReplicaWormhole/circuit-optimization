"""Independent run395 preflight and complete endpoint gate-list audit.

Uses only deterministic NumPy/SciPy matrix and assignment operations. Does not
run the fit, optimize, differentiate, or perturb any circuit parameters.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
FIX_PATH = HERE / "objective_fixtures.py"
FIX_SPEC = importlib.util.spec_from_file_location("label_objective_fixtures", FIX_PATH)
fixtures = importlib.util.module_from_spec(FIX_SPEC)
FIX_SPEC.loader.exec_module(fixtures)
base = fixtures.base
SLOTS = fixtures.LABELS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_family(config: dict):
    spec = importlib.util.spec_from_file_location("full98_label_family_audit", ROOT / config["family_module"])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def assert_static_plan(run_dir: Path, cfg: dict) -> dict:
    script = run_dir / "fit.py"
    if digest(script) != cfg["code_sha256"]:
        raise SystemExit("fit.py does not match the frozen code hash")
    for rel, expected in cfg["dependencies"].items():
        if digest(ROOT / rel) != expected:
            raise SystemExit(f"frozen dependency mismatch: {rel}")
    if scipy.__version__ != cfg["scipy_version"]:
        raise SystemExit(f"SciPy version {scipy.__version__} differs from frozen {cfg['scipy_version']}")
    expected_caps = {"maxiter": 1000, "maxfun": 20000, "maxls": 40,
                     "internal_seconds": 110, "external_seconds": 120, "threads": 1,
                     "nearhit_loss": 1e-18}
    for key, expected in expected_caps.items():
        if cfg.get(key) != expected:
            raise SystemExit(f"frozen fit cap {key} differs from assignment")
    if cfg.get("seed") != 9261620 or cfg.get("initialization") != "exact best98 vector from run394 row seed9261620; no random draw or perturbation":
        raise SystemExit("warm-start source differs from run394 seed9261620")
    parent_path = ROOT / cfg["parent_result"]
    if digest(parent_path) != cfg["dependencies"][cfg["parent_result"]]:
        raise SystemExit("run394 source result hash changed")
    if (run_dir / "result.json").exists():
        raise SystemExit("run395 already has results; preflight-only path refused")
    source = script.read_text()
    code_fragments = (
        "SLOTS = np.array([1] * 6 + [-1] * 4 + [1j] * 3 + [-1j] * 3",
        "cost = np.sum(np.abs(uv[:, None, :] - SLOTS[None, :, None] * u[:, None, :]) ** 2, axis=2)",
        "rows, columns = linear_sum_assignment(cost)",
        "d = torch.tensor(labels, dtype=torch.complex128)",
        "residual = u @ v - d[:, None] * u",
        "gradient = torch.autograd.grad(value, x)[0]",
    )
    for fragment in code_fragments:
        if fragment not in source:
            raise SystemExit(f"frozen implementation lacks required expression: {fragment}")
    call = source[source.index("opt = minimize("):source.index("except BudgetStop:")]
    if "bounds=" in call:
        raise SystemExit("unexpected bounds or coordinate mask in L-BFGS-B call")
    return {
        "run_id": int(run_dir.name), "board_hypothesis_id": 82,
        "fit_script_sha256": digest(script), "config_sha256": digest(run_dir / "config.json"),
        "parent_result_sha256": digest(parent_path), "scipy_version": scipy.__version__,
        "caps": expected_caps,
        "all_98_coordinates_free": "total 98; all coordinates free" in cfg["parameter_order"],
        "static_assignment_review": "16x16 unperturbed raw row costs; repeated exact slots; SciPy linear_sum_assignment; selected D is a constant torch tensor while residual is differentiated",
        "scope": "static hash/input review plus finite fixtures; no optimization, gradients, derivatives, or search",
    }


def unique_label_assignment(costs: np.ndarray, labels: np.ndarray,
                            row_ind: np.ndarray, col_ind: np.ndarray,
                            best_total: float) -> dict:
    baseline = np.empty(16, dtype=complex)
    baseline[row_ind] = labels[col_ind]
    alternatives = []
    for row, col in zip(row_ind, col_ind):
        changed = costs.copy()
        changed[row, col] = np.finfo(float).max
        alt_rows, alt_cols = linear_sum_assignment(changed)
        alt_total = float(changed[alt_rows, alt_cols].sum())
        alt_labels = np.empty(16, dtype=complex)
        alt_labels[alt_rows] = labels[alt_cols]
        if np.array_equal(alt_labels, baseline):
            continue
        alternatives.append(alt_total - best_total)
    scale = max(1.0, abs(best_total))
    tie_tol = 64 * np.finfo(float).eps * scale
    tied = [delta for delta in alternatives if abs(delta) <= tie_tol]
    return {"unique_label_vector_at_input": not bool(tied),
            "distinct_label_alternatives_tested": len(alternatives),
            "tie_tolerance_raw_total_cost": tie_tol,
            "minimum_alternative_cost_gap": min(alternatives) if alternatives else None,
            "tied_distinct_label_alternatives": len(tied)}


def preflight(run_dir: Path, cfg: dict) -> dict:
    report = assert_static_plan(run_dir, cfg)
    parent = json.loads((ROOT / cfg["parent_result"]).read_text())
    row = next(r for r in parent["rows"] if r["seed"] == cfg["seed"])
    x = np.asarray(row["best"], dtype=float)
    if x.shape != (98,) or not np.isfinite(x).all():
        raise SystemExit("parent start must be a finite 98-vector")
    expected_x = np.asarray(json.loads((ROOT / "collaboration/work/verifier98continue/POSTFIT_AUDIT.json").read_text())["starts"]
                            [cfg["seed"] - 9261600]["seed"], dtype=float) if False else x
    # Independent family matrix and complete source-circuit check.
    family = import_family(cfg)
    u = base.coordinate_matrix(x)
    parent_native = ROOT / f"experiments/runs/394/endpoint_{cfg['seed']}_native.json"
    parent_compiled = ROOT / f"experiments/runs/394/endpoint_{cfg['seed']}_compiled.json"
    un, native_gates = base.gate_list(parent_native)
    uc, compiled_gates = base.gate_list(parent_compiled)
    if max(base.phase_error(u, un), base.phase_error(un, uc)) > 2e-12:
        raise SystemExit("run394 warmstart vector/native/compiled matrix mismatch")
    generated_native = family.serialize(x)
    generated_compiled = family.helper.compile_native(generated_native)
    if max(base.phase_error(u, base.gate_list(generated_native)[0]),
           base.phase_error(base.gate_list(generated_native)[0], base.gate_list(generated_compiled)[0])) > 2e-12:
        raise SystemExit("regenerated warmstart native/compiled serialization mismatch")
    if (sum(g["gate"] == "cx" for g in native_gates),
            sum(g["gate"] == "xx_yy" for g in native_gates),
            sum(g["gate"] == "cx" for g in compiled_gates)) != (4, 4, 12):
        raise SystemExit("warmstart circuit has unexpected gate costs")

    v = base.right_shift()
    labels, costs, pairs = fixtures.assign_labels(u, v)
    assigned_cost = float(costs[pairs[0], pairs[1]].sum() / 16)
    direct = u @ v - np.diag(labels) @ u
    b = u @ v @ u.conj().T
    off = b - np.diag(np.diag(b))
    direct_loss = float(np.linalg.norm(direct, "fro") ** 2 / 16)
    conjugated_loss = float(np.linalg.norm(b - np.diag(labels), "fro") ** 2 / 16)
    original_loss = float(np.linalg.norm(off, "fro") ** 2 / 16)
    # Re-run the exact frozen SciPy assignment on an identical cost input.
    r2, c2 = linear_sum_assignment(costs.copy())
    if not (np.array_equal(r2, pairs[0]) and np.array_equal(c2, pairs[1])):
        raise SystemExit("identical frozen start cost input did not repeat its assignment")
    if max(abs(direct_loss-assigned_cost), abs(direct_loss-conjugated_loss)) > 2e-12:
        raise SystemExit("raw Hungarian cost, direct UV-DU residual, and conjugated residual disagree")
    roots, multiplicities = np.unique(labels, return_counts=True)
    counts = {complex(z): int(n) for z, n in zip(roots, multiplicities)}
    if counts != {1+0j: 6, -1+0j: 4, 1j: 3, -1j: 3}:
        raise SystemExit(f"selected D has wrong spectrum multiplicities: {counts}")
    assignment_unique = unique_label_assignment(costs, SLOTS, pairs[0], pairs[1], float(costs[pairs[0], pairs[1]].sum()))
    report.update({
        "warmstart_seed": cfg["seed"], "warmstart_vector_sha256_f64le": hashlib.sha256(np.asarray(x, dtype="<f8").tobytes()).hexdigest(),
        "warmstart_native_sha256": digest(parent_native), "warmstart_compiled_sha256": digest(parent_compiled),
        "warmstart_coordinate_vs_parent_native_error_up_to_phase": base.phase_error(u, un),
        "warmstart_parent_native_vs_compiled_error_up_to_phase": base.phase_error(un, uc),
        "warmstart_native_cx_count": sum(g["gate"] == "cx" for g in native_gates),
        "warmstart_native_xx_yy_count": sum(g["gate"] == "xx_yy" for g in native_gates),
        "warmstart_compiled_cx_count": sum(g["gate"] == "cx" for g in compiled_gates),
        "scipy_assignment_columns": pairs[1].tolist(),
        "selected_label_counts": {str(k): val for k, val in counts.items()},
        "selected_label_vector_unique": assignment_unique,
        "same_input_assignment_repeatable": True,
        "raw_cost_sha256_f64le": hashlib.sha256(np.asarray(costs, dtype="<f8").tobytes()).hexdigest(),
        "assigned_cost_normalized": assigned_cost,
        "direct_label_loss": direct_loss,
        "conjugated_label_loss": conjugated_loss,
        "original_offdiagonal_loss": original_loss,
        "original_offdiagonal_max": float(np.max(np.abs(off))),
        "direct_vs_assigned_abs_difference": abs(direct_loss-assigned_cost),
        "direct_vs_conjugated_abs_difference": abs(direct_loss-conjugated_loss),
        "scope": "one frozen starting point and finite objective checks; no optimizer or derivative",
    })
    return report


def audit_postfit(run_dir: Path, cfg: dict, pre: dict) -> dict:
    result_path = run_dir / "result.json"
    result = json.loads(result_path.read_text())
    if result.get("run_id") != int(run_dir.name) or result.get("config") != cfg:
        raise SystemExit("run ID/config in result differs from frozen reservation")
    if result.get("seed") != cfg["seed"] or not np.array_equal(np.asarray(result["initial"]),
                                                              np.asarray(json.loads((ROOT / cfg["parent_result"]).read_text())["rows"][cfg["seed"]-9261600]["best"])):
        raise SystemExit("result initial point differs from frozen run394 warmstart")
    endpoints = result.get("endpoints", {})
    required = {"initial", "best", "best_original", "last_evaluated"}
    if result.get("optimizer_message") is not None:
        required.add("optimizer_return")
    if not required.issubset(endpoints):
        raise SystemExit(f"missing endpoint records: {required - set(endpoints)}")
    rows = []
    expected_files = set()
    for name, record in endpoints.items():
        point = np.asarray(record["point"], dtype=float)
        if point.shape != (98,) or not np.isfinite(point).all():
            raise SystemExit(f"endpoint {name} point is not a finite 98-vector")
        u = base.coordinate_matrix(point)
        native_path, compiled_path = run_dir / record["native"], run_dir / record["compiled"]
        expected_files.update((native_path.name, compiled_path.name))
        if digest(native_path) != record["native_sha256"] or digest(compiled_path) != record["compiled_sha256"]:
            raise SystemExit(f"endpoint {name} gate-list file hash mismatch")
        un, native_gates = base.gate_list(native_path)
        uc, compiled_gates = base.gate_list(compiled_path)
        coord_err, compile_err = base.phase_error(u, un), base.phase_error(un, uc)
        unitarity_n = float(np.max(np.abs(un.conj().T @ un - np.eye(16))))
        unitarity_c = float(np.max(np.abs(uc.conj().T @ uc - np.eye(16))))
        if max(coord_err, compile_err, unitarity_n, unitarity_c) > 2e-12:
            raise SystemExit(f"endpoint {name} matrix or unitarity check failed")
        counts = (sum(g["gate"] == "cx" for g in native_gates),
                  sum(g["gate"] == "xx_yy" for g in native_gates),
                  sum(g["gate"] == "cx" for g in compiled_gates))
        if counts != (4, 4, 12):
            raise SystemExit(f"endpoint {name} gate costs differ: {counts}")
        pieces = fixtures.objective_components(u)
        columns = record.get("columns")
        if columns is not None and columns != pieces["assignment_columns"]:
            raise SystemExit(f"endpoint {name} saved assignment columns differ under pinned solver")
        label_loss = pieces["label_residual_direct"]
        original = pieces["original_offdiagonal_loss"]
        if "label_loss" in record and abs(float(record["label_loss"])-label_loss)>2e-12:
            raise SystemExit(f"endpoint {name} label residual differs from saved objective")
        if "original_loss" in record and abs(float(record["original_loss"])-original)>2e-12:
            raise SystemExit(f"endpoint {name} original offdiagonal loss differs")
        rows.append({
            "endpoint": name, "label_loss_independent": label_loss,
            "original_offdiagonal_loss_independent": original,
            "original_offdiagonal_max": pieces["original_offdiagonal_max"],
            "assigned_cost_loss": pieces["raw_assignment_cost_normalized"],
            "assignment_columns": pieces["assignment_columns"],
            "assignment_matches_saved": columns is None or columns == pieces["assignment_columns"],
            "coordinate_vs_native_error_up_to_phase": coord_err,
            "native_vs_compiled_error_up_to_phase": compile_err,
            "native_unitarity_max_error": unitarity_n,
            "compiled_unitarity_max_error": unitarity_c,
            "native_cx_count": counts[0], "native_xx_yy_count": counts[1],
            "compiled_cx_count": counts[2],
            "native_sha256": digest(native_path), "compiled_sha256": digest(compiled_path),
            "label_error_identity": pieces["identity_abs_difference"],
            "raw_assignment_cost_error": pieces["cost_abs_difference"],
        })
    endpoint_files = {p.name for p in run_dir.glob("*_native.json")} | {p.name for p in run_dir.glob("*_compiled.json")}
    if endpoint_files != expected_files:
        raise SystemExit(f"unmatched saved endpoint lists: expected {sorted(expected_files)}, got {sorted(endpoint_files)}")
    history = result.get("history", [])
    if len(history) != result.get("evaluations"):
        raise SystemExit("history length differs from evaluation count")
    previous_columns = None
    assignment_changes = 0
    for index, entry in enumerate(history, start=1):
        if entry.get("evaluation") != index or len(entry.get("columns", [])) != 16:
            raise SystemExit("malformed saved assignment history")
        if previous_columns is not None and entry["columns"] != previous_columns:
            assignment_changes += 1
        previous_columns = entry["columns"]
    best_record = endpoints["best"]
    max_original = max(r["original_offdiagonal_max"] for r in rows)
    out = dict(pre)
    out.update({
        "result_sha256": digest(result_path), "termination": result["termination"],
        "elapsed_seconds": result["elapsed_seconds"], "evaluations": result["evaluations"],
        "iterations": result["nit"], "saved_assignment_changes": assignment_changes,
        "best_label_loss": best_record["label_loss"],
        "best_original_loss": best_record["original_loss"],
        "best_original_maxoff": next(r["original_offdiagonal_max"] for r in rows if r["endpoint"] == "best"),
        "endpoints": rows,
        "scope": "independent finite matrix, dynamic-assignment, and original residual audit of every saved endpoint; no exact certificate",
    })
    return out


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: python3 audit_label_run.py RESERVED_RUN_DIR [postfit]")
    run_dir = Path(sys.argv[1]).resolve()
    cfg = json.loads((run_dir / "config.json").read_text())
    if len(sys.argv) == 2:
        report = preflight(run_dir, cfg)
        out = HERE / "PREFIT_APPROVAL.json"
    elif sys.argv[2] == "postfit":
        pre = json.loads((HERE / "PREFIT_APPROVAL.json").read_text())
        report = audit_postfit(run_dir, cfg, pre)
        out = HERE / "POSTFIT_AUDIT.json"
    else:
        raise SystemExit("phase must be 'postfit'")
    report["auditor_script_sha256"] = digest(Path(__file__))
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in report if k != "endpoints"}, indent=2))


if __name__ == "__main__":
    main()
