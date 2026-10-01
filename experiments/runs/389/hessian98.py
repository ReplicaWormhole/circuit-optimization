"""One frozen full-coupled 98-coordinate AD Hessian diagnostic; no fits."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from delete14_search import axis_rotation, axis_to_u3, kron_gate  # noqa: E402
from check_circuit import evaluate, right_shift  # noqa: E402
from search13_adaptive_topology import matrix_for_schedule  # noqa: E402
from native_gate_check import evaluate as native_evaluate  # noqa: E402

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "config.json").read_text())
HELPER_PATH = ROOT / "collaboration/work/verifier_label12/parity12/native_helpers.py"
_spec = importlib.util.spec_from_file_location("native_helpers_98", HELPER_PATH)
helper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(helper)

PAIRS = [(0, 2), (1, 3), (0, 1), (2, 3)]
CX = matrix_for_schedule([(0, 1), (2, 1), (3, 1), (1, 2)])
CX_PAIRS = [(0, 1), (2, 1), (3, 1), (1, 2)]


def matrix(z):
    """Chronology L0 CX01 A1 CX21 B1 CX31 L1 F02 ... L6."""
    local = z[:84].reshape(7, 4, 3)
    ab = z[84:90].reshape(2, 3)
    f = z[90:].reshape(4, 2)
    u = torch.eye(16, dtype=torch.complex128)
    for q in range(4):
        u = kron_gate(axis_rotation(*local[0, q]), q) @ u
    u = CX[0] @ u
    u = kron_gate(axis_rotation(*ab[0]), 1) @ u
    u = CX[1] @ u
    u = kron_gate(axis_rotation(*ab[1]), 1) @ u
    u = CX[2] @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[1, q]), q) @ u
    u = helper.f_matrix(f[0, 0], f[0, 1], *PAIRS[0]) @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[2, q]), q) @ u
    u = helper.f_matrix(f[1, 0], f[1, 1], *PAIRS[1]) @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[3, q]), q) @ u
    u = CX[3] @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[4, q]), q) @ u
    u = helper.f_matrix(f[2, 0], f[2, 1], *PAIRS[2]) @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[5, q]), q) @ u
    u = helper.f_matrix(f[3, 0], f[3, 1], *PAIRS[3]) @ u
    for q in range(4):
        u = kron_gate(axis_rotation(*local[6, q]), q) @ u
    return u


def objective(z):
    u = matrix(z)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    b = u @ v @ u.conj().T
    off = b - torch.diag(torch.diagonal(b))
    return (off.abs() ** 2).sum().real / 16


def serialize(z):
    local = z[:84].reshape(7, 4, 3)
    ab = z[84:90].reshape(2, 3)
    f = z[90:].reshape(4, 2)
    gates = []
    def add_layer(k):
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q, **axis_to_u3(*local[k, q])})
    add_layer(0)
    gates.append({"gate": "cx", "control": 0, "target": 1})
    gates.append({"gate": "u3", "qubit": 1, **axis_to_u3(*ab[0])})
    gates.append({"gate": "cx", "control": 2, "target": 1})
    gates.append({"gate": "u3", "qubit": 1, **axis_to_u3(*ab[1])})
    gates.append({"gate": "cx", "control": 3, "target": 1})
    add_layer(1)
    gates.append(helper.native_gate(*PAIRS[0], *f[0]))
    add_layer(2)
    gates.append(helper.native_gate(*PAIRS[1], *f[1]))
    add_layer(3)
    gates.append({"gate": "cx", "control": 1, "target": 2})
    add_layer(4)
    gates.append(helper.native_gate(*PAIRS[2], *f[2]))
    add_layer(5)
    gates.append(helper.native_gate(*PAIRS[3], *f[3]))
    add_layer(6)
    return {"n": 4, "gates": gates}


def remap116_to98(x116):
    """Boundary-absorb A/B controls; retain A1/B1; return frozen 98 point."""
    old = np.asarray(x116[:108], dtype=float).reshape(9, 4, 3)
    def su2(r):
        x, y, z = r
        length = np.linalg.norm(r)
        c = math.cos(length / 2)
        s = 0.5 * np.sinc(length / (2 * math.pi))
        return np.array([[c - 1j*z*s, (-1j*x-y)*s], [(-1j*x+y)*s, c + 1j*z*s]])
    def rotvec(m):
        c = float(np.clip(np.trace(m).real / 2, -1.0, 1.0))
        v = np.array([
            -(m[0, 1].imag + m[1, 0].imag) / 2,
            (m[1, 0].real - m[0, 1].real) / 2,
            (m[1, 1].imag - m[0, 0].imag) / 2,
        ])
        s = float(np.linalg.norm(v))
        if s < 1e-12:
            if c > 0:
                return np.zeros(3)
            raise ValueError("unstable pi-rotation remap")
        return (2 * math.atan2(s, c) / s) * v
    mapped = np.asarray([old[0].copy(), old[3].copy(), *[a.copy() for a in old[4:9]]])
    # Product order follows the chronological absorption identities in ROUND_13.
    mapped[0, 2] = rotvec(su2(old[1, 2]) @ su2(old[0, 2]))
    mapped[0, 3] = rotvec(su2(old[2, 3]) @ su2(old[1, 3]) @ su2(old[0, 3]))
    mapped[1, 0] = rotvec(su2(old[3, 0]) @ su2(old[2, 0]) @ su2(old[1, 0]))
    mapped[1, 2] = rotvec(su2(old[3, 2]) @ su2(old[2, 2]))
    locals98 = np.asarray(mapped).reshape(84)
    ab = np.concatenate((old[1, 1], old[2, 1]))
    return np.concatenate((locals98, ab, np.asarray(x116[108:], dtype=float)))


def main():
    if (HERE / "result.json").exists():
        raise SystemExit("guard: result already exists; refusing a repeated run")
    for rel, digest in CFG["dependencies"].items():
        actual = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        if actual != digest:
            raise SystemExit(f"dependency hash changed: {rel}")
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != CFG["code_sha256"]:
        raise SystemExit("run script hash differs from frozen configuration")
    for rel, digest in CFG["source_hashes"].items():
        source = HERE / rel if rel == "run387_frozen_parent.json" else ROOT / rel
        if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            raise SystemExit(f"frozen source hash changed: {rel}")
    point = json.loads((HERE / "run387_frozen_parent.json").read_text())["config"]["point116"]
    point98 = remap116_to98(point)
    if point98.shape != (98,):
        raise SystemExit(f"bad embedded point shape {point98.shape}")
    torch.set_num_threads(1)
    x = torch.tensor(point98, dtype=torch.float64, requires_grad=True)
    v = np.asarray(right_shift(4), dtype=np.complex128)
    t0 = time.monotonic()
    with threadpool_limits(limits=1):
        loss = objective(x)
        grad = torch.autograd.grad(loss, x, create_graph=False)[0]
        hessian = torch.autograd.functional.hessian(objective, x, vectorize=False)
    elapsed = time.monotonic() - t0
    h = hessian.detach().cpu().numpy()
    hsym = (h + h.T) / 2
    eigvals, eigvecs = np.linalg.eigh(hsym)
    for j in range(eigvecs.shape[1]):
        i = int(np.argmax(np.abs(eigvecs[:, j])))
        if eigvecs[i, j] < 0:
            eigvecs[:, j] *= -1
    u = matrix(x).detach().cpu().numpy()
    native = serialize(point98)
    compiled = helper.compile_native(native)
    native_path, compiled_path = HERE / "base_native.json", HERE / "base_compiled.json"
    native_path.write_text(json.dumps(native, indent=2) + "\n")
    compiled_path.write_text(json.dumps(compiled, indent=2) + "\n")
    ncheck = native_evaluate(native)
    ccheck = evaluate(compiled)
    out = {
        "run_id": CFG.get("run_id") or (int(HERE.name) if HERE.name.isdigit() else None),
        "seed": CFG["seed"], "parent_run": 385,
        "embedding_parent_run": 387, "elapsed_seconds": elapsed,
        "threads": 1, "loss": float(loss.detach()), "gradient": grad.detach().cpu().numpy().tolist(),
        "gradient_norm": float(torch.linalg.vector_norm(grad)),
        "hessian": h.tolist(), "hessian_symmetry_max": float(np.max(np.abs(h - h.T))),
        "eigenvalues": eigvals.tolist(), "eigenvectors_columns": eigvecs.tolist(),
        "least_eigenvalue": float(eigvals[0]),
        "least_eigenpair_residual": float(np.linalg.norm(hsym @ eigvecs[:, 0] - eigvals[0] * eigvecs[:, 0])),
        "negative_trigger_threshold": CFG["negative_curvature_trigger"],
        "trigger_fires": bool(eigvals[0] < CFG["negative_curvature_trigger"]),
        "point98": point98.tolist(), "native_path": native_path.name,
        "compiled_path": compiled_path.name,
        "native_sha256": hashlib.sha256(native_path.read_bytes()).hexdigest(),
        "compiled_sha256": hashlib.sha256(compiled_path.read_bytes()).hexdigest(),
        "native_checker": ncheck, "compiled_checker": ccheck,
        "direct_matrix_unitarity_max": float(np.max(np.abs(u.conj().T @ u - np.eye(16)))),
        "scope": "one floating-point full98 Hessian diagnostic; no optimizer, no fits, no exclusion or certificate",
    }
    (HERE / "result.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("run_id", "elapsed_seconds", "loss", "gradient_norm", "least_eigenvalue", "least_eigenpair_residual", "trigger_fires", "native_checker", "compiled_checker")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
