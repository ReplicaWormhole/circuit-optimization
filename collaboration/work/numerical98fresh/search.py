"""One frozen, single-start full98 L-BFGS-B fit; coordinator runs after reservation."""
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import torch
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
CFG = json.loads((HERE / "config.json").read_text())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_family():
    path = ROOT / CFG["family_module"]
    spec = importlib.util.spec_from_file_location("frozen_full98_family", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    try:
        run_id = int(HERE.name)
    except ValueError as error:
        raise SystemExit("guard: run search from its reserved numeric workspace") from error
    result_path = HERE / "result.json"
    if result_path.exists():
        raise SystemExit("guard: result already exists; refusing a repeated run")
    if sha256(__file__) != CFG["code_sha256"]:
        raise SystemExit("guard: search.py differs from its frozen code hash")
    for rel, expected in CFG["source_hashes"].items():
        path = HERE / rel if rel == "config.json" else ROOT / rel
        if sha256(path) != expected:
            raise SystemExit(f"guard: frozen source changed: {rel}")
    family = load_family()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    x0 = np.asarray(CFG["initial_vector"], dtype=np.float64)
    if x0.shape != (98,) or not np.isfinite(x0).all():
        raise SystemExit("guard: frozen initial vector must contain 98 finite values")
    if np.linalg.norm(x0[84:87]) <= 1e-12 or np.linalg.norm(x0[87:90]) <= 1e-12:
        raise SystemExit("guard: A1 and B1 must both be nonzero")
    # Check the seed rule against the literal frozen vector before any fit.
    rng = np.random.default_rng(CFG["seed"])
    reproduced = rng.normal(0.0, 0.3, size=98)
    if not np.array_equal(reproduced, x0):
        raise SystemExit("guard: frozen vector does not reproduce from the stated RNG rule")

    t0 = time.monotonic()
    deadline = t0 + CFG["internal_timeout_seconds"]
    history = []
    best = {"loss": float("inf"), "vector": None, "gradient": None}
    eval_count = 0
    stop_reason = "optimizer_returned"
    timeout_error = None

    def value_and_gradient(vector):
        nonlocal eval_count, best
        if time.monotonic() >= deadline:
            raise TimeoutError("internal monotonic 110-second fit deadline reached")
        if eval_count >= CFG["maxfun"]:
            raise EvaluationBudget("frozen maxfun=10000 evaluation budget reached")
        point = np.asarray(vector, dtype=np.float64)
        tx = torch.tensor(point, dtype=torch.float64, requires_grad=True)
        loss = family.objective(tx)
        gradient = torch.autograd.grad(loss, tx)[0]
        loss_value = float(loss.detach())
        gradient_value = gradient.detach().cpu().numpy().copy()
        eval_count += 1
        elapsed = time.monotonic() - t0
        grad_norm = float(np.linalg.norm(gradient_value))
        history.append({"evaluation": eval_count, "elapsed_seconds": elapsed,
                        "loss": loss_value, "gradient_norm": grad_norm})
        if loss_value < best["loss"]:
            best = {"loss": loss_value, "vector": point.copy(),
                    "gradient": gradient_value.copy(), "evaluation": eval_count,
                    "elapsed_seconds": elapsed}
        if time.monotonic() >= deadline:
            raise TimeoutError("internal monotonic 110-second fit deadline reached after saving evaluation")
        return loss_value, gradient_value

    optimizer_result = None
    started_at = time.time()
    with threadpool_limits(limits=1):
        try:
            optimizer_result = minimize(
                value_and_gradient, x0.copy(), jac=True, method="L-BFGS-B",
                options={"maxiter": CFG["maxiter"], "maxfun": CFG["maxfun"],
                         "maxls": CFG["maxls"], "ftol": CFG["ftol"],
                         "gtol": CFG["gtol"]})
        except TimeoutError as error:
            stop_reason = "internal_timeout"
            timeout_error = str(error)
        except EvaluationBudget as error:
            stop_reason = "maxfun_budget"
            timeout_error = str(error)
    if best["vector"] is None:
        raise RuntimeError("fit ended without an evaluated point")

    initial_native = family.serialize(x0)
    initial_compiled = family.helper.compile_native(initial_native)
    best_native = family.serialize(best["vector"])
    best_compiled = family.helper.compile_native(best_native)
    gate_lists = {
        "initial_native.json": initial_native,
        "initial_compiled.json": initial_compiled,
        "best_native.json": best_native,
        "best_compiled.json": best_compiled,
    }
    for name, candidate in gate_lists.items():
        (HERE / name).write_text(json.dumps(candidate, indent=2) + "\n")
    native_check = family.native_evaluate(best_native)
    compiled_check = family.evaluate(best_compiled)
    initial_native_check = family.native_evaluate(initial_native)
    initial_compiled_check = family.evaluate(initial_compiled)
    elapsed_total = time.monotonic() - t0
    versions = {
        "python": platform.python_version(), "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "torch": torch.__version__,
        "threadpoolctl": __import__("threadpoolctl").__version__}
    out = {
        "hypothesis_id": CFG["hypothesis_id"], "run_id": run_id,
        "seed": CFG["seed"], "initialization_rule": CFG["initialization_rule"],
        "parameter_order": CFG["parameter_order"],
        "coordinates_free": CFG["parameters"], "A1_vector": x0[84:87].tolist(),
        "B1_vector": x0[87:90].tolist(),
        "native_schedule": CFG["native_schedule"],
        "native_entanglers": CFG["native_entanglers"],
        "compiled_cnot_count": CFG["compiled_cnot_count"],
        "optimizer": {
            "method": "L-BFGS-B", "jacobian": "PyTorch autograd",
            "options": {k: CFG[k] for k in ("maxiter", "maxfun", "maxls", "ftol", "gtol")},
            "stop_reason": stop_reason, "timeout_error": timeout_error,
            "success": None if optimizer_result is None else bool(optimizer_result.success),
            "message": None if optimizer_result is None else str(optimizer_result.message),
            "nit": None if optimizer_result is None else int(optimizer_result.nit),
            "nfev": None if optimizer_result is None else int(optimizer_result.nfev),
            "njev": None if optimizer_result is None else int(optimizer_result.njev),
            "reported_fun": None if optimizer_result is None else float(optimizer_result.fun),
            "reported_x": None if optimizer_result is None else np.asarray(optimizer_result.x).tolist()},
        "stop_reason": stop_reason,
        "counts": {"objective_gradient_evaluations": eval_count,
                   "starts": 1, "restarts": 0, "threads": 1},
        "initial_loss": history[0]["loss"],
        "best_loss": best["loss"], "best_evaluation": best["evaluation"],
        "best_gradient_norm": float(np.linalg.norm(best["gradient"])),
        "best_gradient": best["gradient"].tolist(),
        "initial_vector": x0.tolist(), "best_vector": best["vector"].tolist(),
        "initial_native": initial_native, "initial_compiled": initial_compiled,
        "best_native": best_native, "best_compiled": best_compiled,
        "initial_native_checker": initial_native_check,
        "initial_compiled_checker": initial_compiled_check,
        "best_native_checker": native_check, "best_compiled_checker": compiled_check,
        "history": history,
        "runtime": {"started_unix": started_at, "elapsed_seconds": elapsed_total,
                    "internal_limit_seconds": CFG["internal_timeout_seconds"],
                    "external_limit_seconds": CFG["external_timeout_seconds"]},
        "versions": versions,
        "source_hashes": CFG["source_hashes"], "code_sha256": CFG["code_sha256"],
        "scope": "one numerical start in the fixed 4-CX+4-XX/YY schedule; coordinate-wise fresh seed only; no exact certification or family/lower-bound exclusion"
    }
    result_path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("run_id", "stop_reason", "counts", "initial_loss", "best_loss",
                       "best_gradient_norm", "runtime")}, indent=2), flush=True)


class EvaluationBudget(Exception):
    """Raised before an objective evaluation would exceed frozen maxfun."""


if __name__ == "__main__":
    main()
