"""One frozen full98 fit for the reordered Bell-decoder word."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import importlib.util
import json
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


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def import_family(path):
    spec = importlib.util.spec_from_file_location("bellprefix_family", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def main():
    cfg = json.loads((HERE / "config.json").read_text())
    marker = HERE / "execution_started.json"
    if (HERE / "result.json").exists() or marker.exists():
        raise RuntimeError("Refusing to repeat or restart this frozen single start")
    if digest(__file__) != cfg["code_sha256"]:
        raise RuntimeError("search.py differs from frozen code hash")
    if digest(HERE / "family.py") != cfg["family_sha256"]:
        raise RuntimeError("family.py differs from frozen family hash")
    for rel, expected in cfg["dependencies"].items():
        if digest(ROOT / rel) != expected:
            raise RuntimeError(f"frozen dependency changed: {rel}")
    if scipy.__version__ != cfg["scipy_version"]:
        raise RuntimeError(f"SciPy version mismatch: {scipy.__version__}")
    marker.write_text(json.dumps({"run_id": int(HERE.name) if HERE.name.isdigit() else None,
                                 "code_sha256": cfg["code_sha256"],
                                 "started_unix": time.time()}, indent=2) + "\n")

    family = import_family(HERE / "family.py")
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    base, noise, initial = family.initial_point(cfg["seed"], cfg["noise_sigma"])
    initial_native = family.serialize(initial)
    initial_compiled = family.compile_native(initial_native)
    write_json(HERE / "initial_native.json", initial_native)
    write_json(HERE / "initial_compiled.json", initial_compiled)

    started = time.monotonic()
    deadline = started + cfg["internal_seconds"]
    calls = 0
    history = []
    best = {"loss": None, "x": initial.copy(), "gradient": None,
            "evaluation": None, "evaluated": False}
    optimizer_result = None
    termination = "optimizer_returned"

    class BudgetStop(Exception):
        pass

    def objective_and_gradient(z):
        nonlocal calls, best
        if time.monotonic() >= deadline:
            raise BudgetStop("internal wall-clock deadline reached")
        if calls >= cfg["maxfun"]:
            raise BudgetStop("hard objective-evaluation cap reached")
        x = np.asarray(z, dtype=np.float64).copy()
        variable = torch.tensor(x, dtype=torch.float64, requires_grad=True)
        loss_tensor = family.objective(variable)
        gradient = torch.autograd.grad(loss_tensor, variable)[0].detach().cpu().numpy()
        loss = float(loss_tensor.detach().cpu())
        if not np.isfinite(loss) or not np.isfinite(gradient).all():
            raise FloatingPointError("nonfinite objective or gradient")
        calls += 1
        history.append({"evaluation": calls, "loss": loss,
                        "gradient_norm": float(np.linalg.norm(gradient))})
        if best["loss"] is None or loss < best["loss"]:
            best = {"loss": loss, "x": x, "gradient": gradient.copy(),
                    "evaluation": calls, "evaluated": True}
        return loss, gradient

    with threadpool_limits(limits=1):
        try:
            optimizer_result = minimize(
                objective_and_gradient, initial, jac=True, method="L-BFGS-B",
                options={"maxiter": cfg["maxiter"], "maxfun": cfg["maxfun"],
                         "maxls": cfg["maxls"], "ftol": cfg["ftol"],
                         "gtol": cfg["gtol"]})
        except BudgetStop as error:
            termination = "budget_stop: " + str(error)
        except FloatingPointError as error:
            termination = "nonfinite_evaluation: " + str(error)

    best_native = family.serialize(best["x"])
    best_compiled = family.compile_native(best_native)
    write_json(HERE / "best_native.json", best_native)
    write_json(HERE / "best_compiled.json", best_compiled)
    with threadpool_limits(limits=1):
        native_checker = family.evaluate_native(best_native)
        compiled_checker = family.evaluate_compiled(best_compiled)
        initial_native_checker = family.evaluate_native(initial_native)
        initial_compiled_checker = family.evaluate_compiled(initial_compiled)
        best_loss_recomputed = None
        final_recomputation_count = 0
        if best["evaluated"]:
            best_loss_recomputed = float(family.objective(
                torch.tensor(best["x"], dtype=torch.float64)).detach())
            final_recomputation_count = 1
    result = {
        "run_id": int(HERE.name) if HERE.name.isdigit() else None,
        "config": cfg,
        "seed": cfg["seed"],
        "initialization": cfg["initialization"],
        "bell_decoder_base98": base.tolist(),
        "gaussian_perturbation98": noise.tolist(),
        "initial98": initial.tolist(),
        "initial_native_sha256": digest(HERE / "initial_native.json"),
        "initial_compiled_sha256": digest(HERE / "initial_compiled.json"),
        "best98": best["x"].tolist(),
        "best_evaluated": best["evaluated"],
        "best_evaluation": best["evaluation"],
        "best_loss": best["loss"],
        "best_loss_recomputed": best_loss_recomputed,
        "final_recomputation_count_outside_optimizer_maxfun": final_recomputation_count,
        "best_gradient": None if best["gradient"] is None else best["gradient"].tolist(),
        "best_gradient_norm": None if best["gradient"] is None else float(np.linalg.norm(best["gradient"])),
        "best_native_sha256": digest(HERE / "best_native.json"),
        "best_compiled_sha256": digest(HERE / "best_compiled.json"),
        "history": history,
        "evaluations": calls,
        "elapsed_seconds": time.monotonic() - started,
        "termination": termination,
        "optimizer": None if optimizer_result is None else {
            "success": bool(optimizer_result.success),
            "message": str(optimizer_result.message),
            "nit": int(optimizer_result.nit),
            "nfev": int(optimizer_result.nfev),
            "njev": int(optimizer_result.njev),
            "fun": float(optimizer_result.fun),
            "jac": np.asarray(optimizer_result.jac).tolist(),
        },
        "initial_native_checker": initial_native_checker,
        "initial_compiled_checker": initial_compiled_checker,
        "best_native_checker": native_checker,
        "best_compiled_checker": compiled_checker,
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__,
                     "torch": torch.__version__},
        "scope": "One bounded fit in the reordered full98 family; not an exact certificate, family exclusion, or lower bound",
    }
    write_json(HERE / "result.json", result)
    print(json.dumps({"run_id": result["run_id"], "evaluations": calls,
                      "best_loss": best["loss"], "termination": termination,
                      "elapsed_seconds": result["elapsed_seconds"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
