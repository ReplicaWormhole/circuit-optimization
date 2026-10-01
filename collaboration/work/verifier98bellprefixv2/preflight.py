"""Three-point independent finite preflight for the reordered full98 family.
Run only from its separately reserved numeric workspace after board linkage and
root authorization. No objective, derivative, optimizer, or search is evaluated.
"""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm
import torch
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def import_family(family_path):
    spec = importlib.util.spec_from_file_location("frozen_bellprefix_family", family_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def su2(axis):
    """Independent exp(-i axis.sigma/2) NumPy formula."""
    x, y, z = map(float, axis)
    theta = math.sqrt(x * x + y * y + z * z)
    if theta == 0.0:
        return I2.copy()
    return (math.cos(theta / 2) * I2
            - 1j * math.sin(theta / 2) / theta * (x * X + y * Y + z * Z))


def embed(one, q):
    factors = [one if i == q else I2 for i in range(4)]
    out = factors[0]
    for factor in factors[1:]:
        out = np.kron(out, factor)
    return out


def cx(control, target):
    out = np.zeros((16, 16), dtype=complex)
    for col in range(16):
        row = col
        if (col >> (3 - control)) & 1:
            row ^= 1 << (3 - target)
        out[row, col] = 1.0
    return out


def f_gate(a, b, q0, q1):
    xx = embed(X, q0) @ embed(X, q1)
    yy = embed(Y, q0) @ embed(Y, q1)
    return expm(1j * (float(a) * xx + float(b) * yy))


def coordinate_matrix(point):
    """Independent NumPy realization of the frozen chronological 98-word."""
    z = np.asarray(point, dtype=float)
    if z.shape != (98,):
        raise ValueError("expected a 98-coordinate point")
    local = z[:84].reshape(7, 4, 3)
    ab = z[84:90].reshape(2, 3)
    f = z[90:].reshape(4, 2)
    u = np.eye(16, dtype=complex)

    def layer(k):
        nonlocal u
        for q in range(4):
            u = embed(su2(local[k, q]), q) @ u

    def apply_cx(c, t):
        nonlocal u
        u = cx(c, t) @ u

    layer(0)
    u = f_gate(*f[0], 0, 2) @ u
    layer(1)
    u = f_gate(*f[1], 1, 3) @ u
    layer(2)
    apply_cx(0, 1)
    u = embed(su2(ab[0]), 1) @ u
    apply_cx(2, 1)
    u = embed(su2(ab[1]), 1) @ u
    apply_cx(3, 1)
    layer(3)
    apply_cx(1, 2)
    layer(4)
    u = f_gate(*f[2], 0, 1) @ u
    layer(5)
    u = f_gate(*f[3], 2, 3) @ u
    layer(6)
    return u


def one_qubit_from_gate(g):
    kind = g["gate"].lower()
    if kind == "u3":
        t, p, l = (float(g[k]) for k in ("theta", "phi", "lam"))
        c, s = math.cos(t / 2), math.sin(t / 2)
        return np.array([[c, -np.exp(1j * l) * s],
                         [np.exp(1j * p) * s, np.exp(1j * (p + l)) * c]], dtype=complex)
    t = float(g["theta"])
    c, s = math.cos(t / 2), math.sin(t / 2)
    if kind == "rx":
        return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)
    if kind == "ry":
        return np.array([[c, -s], [s, c]], dtype=complex)
    if kind == "rz":
        return np.diag([np.exp(-0.5j * t), np.exp(0.5j * t)])
    raise ValueError(f"unsupported serialized one-qubit gate {kind}")


def gate_matrix(g):
    kind = g["gate"].lower()
    if kind == "cx":
        return cx(int(g["control"]), int(g["target"]))
    if kind == "xx_yy":
        q0, q1 = map(int, g["qubits"])
        return f_gate(float(g["a"]), float(g["b"]), q0, q1)
    return embed(one_qubit_from_gate(g), int(g["qubit"]))


def gate_list_matrix(payload):
    if payload.get("n") != 4 or not isinstance(payload.get("gates"), list):
        raise ValueError("serialized circuit must declare four wires and a gates list")
    u = np.eye(16, dtype=complex)
    for gate in payload["gates"]:
        u = gate_matrix(gate) @ u
    return u


def phase_error(a, b):
    overlap = np.vdot(b, a)
    if abs(overlap) == 0.0:
        return float("inf")
    return float(np.max(np.abs(a - overlap / abs(overlap) * b)))


def decoder_pair_expected():
    """Direct pair-Bell decoder E_pair^dagger in q0,q1,q2,q3 order."""
    bell = np.array([[1, 0, 1, 0], [0, 1, 0, 1],
                     [0, 1, 0, -1], [1, 0, -1, 0]], dtype=complex) / math.sqrt(2)
    regroup = np.zeros((16, 16), dtype=complex)
    for standard in range(16):
        bits = [(standard >> (3 - q)) & 1 for q in range(4)]
        grouped = (((bits[0] << 1) | bits[2]) << 2) | ((bits[1] << 1) | bits[3])
        regroup[grouped, standard] = 1.0
    return regroup.T @ np.kron(bell.conj().T, bell.conj().T) @ regroup


def expected_base():
    x = np.zeros(98, dtype=np.float64)
    local = x[:84].reshape(7, 4, 3)
    f = x[90:].reshape(4, 2)
    dec = 2 * math.pi / (3 * math.sqrt(3)) * np.array([-1.0, 1.0, -1.0])
    local[0, 0] = dec
    local[0, 1] = dec
    local[0, 2] = [-math.pi / 2, 0.0, 0.0]
    local[0, 3] = [-math.pi / 2, 0.0, 0.0]
    f[0] = [-math.pi / 4, 0.0]
    f[1] = [-math.pi / 4, 0.0]
    local[1, 0] = [math.pi, 0.0, 0.0]
    local[2, 1] = [math.pi, 0.0, 0.0]
    return x


def validate_serialized_word(candidate, numeric_config):
    gates = candidate["gates"]
    cursor = 0

    def take_layer(name):
        nonlocal cursor
        for q in range(4):
            if cursor >= len(gates) or gates[cursor].get("gate") != "u3" or gates[cursor].get("qubit") != q:
                raise SystemExit(f"serialized local layer {name} is not chronological q0..q3")
            cursor += 1

    def take_entangler(kind, first, second):
        nonlocal cursor
        if cursor >= len(gates):
            raise SystemExit("serialized gate list ended before its declared word")
        gate = gates[cursor]
        if kind == "cx":
            ok = gate.get("gate") == "cx" and gate.get("control") == first and gate.get("target") == second
        else:
            ok = gate.get("gate") == "xx_yy" and gate.get("qubits") == [first, second]
        if not ok:
            raise SystemExit(f"serialized {kind} order mismatch at gate {cursor}")
        cursor += 1

    def take_special(name, qubit):
        nonlocal cursor
        if cursor >= len(gates) or gates[cursor].get("gate") != "u3" or gates[cursor].get("qubit") != qubit:
            raise SystemExit(f"serialized {name} local rotation is misplaced")
        cursor += 1

    for token in numeric_config["native_word"]:
        kind = token[0]
        if kind == "layer":
            take_layer(token[1])
        elif kind in ("cx", "xx_yy"):
            take_entangler(kind, token[1], token[2])
        elif kind in ("A1", "B1"):
            take_special(kind, token[1])
        else:
            raise SystemExit(f"unknown frozen word token {token}")
    if cursor != len(gates):
        raise SystemExit(f"serialized gate list has {len(gates)-cursor} undeclared extra gates")
    if sum(g["gate"] == "cx" for g in gates) != 4 or sum(g["gate"] == "xx_yy" for g in gates) != 4:
        raise SystemExit("serialized native counts are not four CX plus four XX/YY")
    if sum(g["gate"] == "u3" for g in gates) != 30:
        raise SystemExit("serialized U3 count is not 28 layer rotations plus A1/B1")
    compiled = frozen_family.compile_native(candidate)
    if sum(g["gate"] == "cx" for g in compiled["gates"]) != 12:
        raise SystemExit("compiled circuit does not contain twelve CX")
    if any(g["gate"] == "xx_yy" for g in compiled["gates"]):
        raise SystemExit("compiled circuit retains native XX/YY gates")
    return compiled


def main():
    config_path = HERE / "preflight_config.json"
    frozen = json.loads(config_path.read_text())
    try:
        run_id = int(HERE.name)
    except ValueError as exc:
        raise SystemExit("preflight checker must run from a reserved numeric run directory") from exc
    if frozen.get("run_id") is not None and frozen.get("run_id") != run_id:
        raise SystemExit("frozen preflight run ID does not match workspace")
    numeric_workspace = frozen.get("numeric_fit_workspace")
    if not isinstance(numeric_workspace, str) or not numeric_workspace.startswith("experiments/runs/"):
        raise SystemExit("preflight freeze must name the exact reserved numerical workspace")
    fit_dir = ROOT / numeric_workspace
    family_path = fit_dir / "family.py"
    search_path = fit_dir / "search.py"
    numeric_config_path = fit_dir / "config.json"
    if (HERE / "preflight_result.json").exists():
        raise SystemExit("refusing to repeat this bounded preflight")
    if sha(Path(__file__)) != frozen["checker_sha256"]:
        raise SystemExit("preflight checker hash mismatch")
    if frozen.get("external_seconds") != 30 or frozen.get("threads") != 1:
        raise SystemExit("preflight runtime cap differs from freeze")
    if (frozen.get("max_coordinate_points") != 3 or frozen.get("full_family_matrix_builds") != 6
            or frozen.get("serialized_gate_list_matrix_builds") != 6
            or frozen.get("extra_prefix_and_direct_decoder_matrix_builds") != 2):
        raise SystemExit("preflight matrix-build accounting differs from freeze")
    for rel, expected_hash in frozen["source_hashes"].items():
        if sha(ROOT / rel) != expected_hash:
            raise SystemExit(f"frozen preflight dependency changed: {rel}")
    numeric_config = json.loads(numeric_config_path.read_text())
    if sha(numeric_config_path) != frozen["numeric_config_sha256"]:
        raise SystemExit("numerical fit config hash mismatch")
    if sha(family_path) != frozen["family_sha256"] or sha(search_path) != frozen["search_sha256"]:
        raise SystemExit("numerical family/search hash mismatch")
    if scipy.__version__ != frozen["scipy_version"] or torch.__version__ != frozen["torch_version"]:
        raise SystemExit("pinned NumPy/SciPy/Torch runtime differs")
    if numeric_config["native_word"] != frozen["native_word"]:
        raise SystemExit("declared numerical word differs from frozen preflight word")
    noise_sigma = numeric_config.get("noise_sigma", numeric_config["initialization"]["noise_sigma"])
    if numeric_config["seed"] != frozen["seed"] or noise_sigma != frozen["noise_sigma"]:
        raise SystemExit("seed/sigma mismatch")
    if "noise_sigma" in numeric_config and noise_sigma != numeric_config["initialization"]["noise_sigma"]:
        raise SystemExit("top-level and nested noise sigma disagree")
    if numeric_config["maxiter"] != 1000 or numeric_config["maxfun"] != 20000 or numeric_config["threads"] != 1:
        raise SystemExit("fit bounds differ from frozen config")
    if numeric_config.get("initialization", {}).get("coordinates_perturbed") != "all98" or numeric_config.get("initialization", {}).get("restarts") != 0:
        raise SystemExit("initializer is not one all-98-coordinate start with no restarts")
    if numeric_config.get("native_entanglers") != 8 or numeric_config.get("compiled_cnot_cost") != 12:
        raise SystemExit("declared native or compiled entangler count differs")
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    global frozen_family
    frozen_family = import_family(family_path)

    # Three deterministic points only: exact Bell prefix, frozen RNG start, and an all-nonzero probe.
    base = expected_base()
    if not np.array_equal(frozen_family.bell_decoder_point(), base):
        raise SystemExit("Bell decoder base coordinates differ from independent formula")
    rng = np.random.default_rng(int(numeric_config["seed"]))
    noise = rng.normal(0.0, float(noise_sigma), size=98)
    initial = base + noise
    family_base, family_noise, family_initial = frozen_family.initial_point(
        numeric_config["seed"], noise_sigma)
    if not (np.array_equal(family_base, base) and np.array_equal(family_noise, noise)
            and np.array_equal(family_initial, initial)):
        raise SystemExit("seeded initial point does not match independent default_rng construction")
    if noise.shape != (98,) or not np.isfinite(initial).all():
        raise SystemExit("invalid seeded 98-coordinate start")
    probe = np.linspace(-0.98, 0.98, 98, dtype=np.float64)
    if np.any(probe == 0) or not np.isfinite(probe).all():
        raise SystemExit("deterministic probe must be finite and all nonzero")

    points = [("bell_decoder_base", base), ("seeded_initial", initial), ("all_nonzero_probe", probe)]
    rows = []
    with threadpool_limits(limits=1), torch.no_grad():
        for name, point in points:
            if point.shape != (98,):
                raise SystemExit(f"bad point shape for {name}")
            torch_matrix = frozen_family.matrix(torch.as_tensor(point, dtype=torch.float64)).detach().cpu().numpy()
            numpy_matrix = coordinate_matrix(point)
            native = frozen_family.serialize(point)
            compiled = validate_serialized_word(native, numeric_config)
            native_matrix = gate_list_matrix(native)
            compiled_matrix = gate_list_matrix(compiled)
            errors = {
                "torch_vs_independent_numpy_up_to_phase": phase_error(torch_matrix, numpy_matrix),
                "native_serialization_vs_numpy_up_to_phase": phase_error(native_matrix, numpy_matrix),
                "compiled_vs_native_up_to_phase": phase_error(compiled_matrix, native_matrix),
                "torch_unitarity_max": float(np.max(np.abs(torch_matrix.conj().T @ torch_matrix - np.eye(16)))),
                "numpy_unitarity_max": float(np.max(np.abs(numpy_matrix.conj().T @ numpy_matrix - np.eye(16)))),
            }
            if max(errors.values()) > frozen["tolerance"]:
                raise SystemExit(f"independent matrix/serializer mismatch at {name}: {errors}")
            row = {"name": name, "point_sha256_f64le": hashlib.sha256(np.asarray(point, dtype="<f8").tobytes()).hexdigest(),
                   "coordinate_point_evaluations": 1, "native_gate_count": len(native["gates"]),
                   "native_cx_count": sum(g["gate"] == "cx" for g in native["gates"]),
                   "native_xx_yy_count": sum(g["gate"] == "xx_yy" for g in native["gates"]),
                   "compiled_cx_count": sum(g["gate"] == "cx" for g in compiled["gates"]),
                   "errors": errors}
            if name == "bell_decoder_base":
                prefix = {"n": 4, "gates": native["gates"][:14]}
                prefix_matrix = gate_list_matrix(prefix)
                direct_decoder = decoder_pair_expected()
                row["decoder_prefix_direct_E_pair_dagger_phase_error"] = phase_error(prefix_matrix, direct_decoder)
                if row["decoder_prefix_direct_E_pair_dagger_phase_error"] > frozen["tolerance"]:
                    raise SystemExit("serialized L0 F02 L1 F13 L2 prefix is not E_pair^dagger")
            rows.append(row)
    result = {
        "run_id": run_id,
        "board_hypothesis_id": frozen["board_hypothesis_id"],
        "schema": "reordered-full98-static-independent-preflight-v1",
        "checker_sha256": sha(Path(__file__)),
        "config_sha256": sha(config_path),
        "source_hashes": frozen["source_hashes"],
        "coordinate_point_evaluation_count": 3,
        "full_family_matrix_builds": 6,
        "serialized_gate_list_matrix_builds": 6,
        "extra_prefix_and_direct_decoder_matrix_builds": 2,
        "points": rows,
        "seed": int(numeric_config["seed"]),
        "sigma": float(noise_sigma),
        "all98_coordinates_free_by_declared_slices": True,
        "native_counts": {"cx": 4, "xx_yy": 4},
        "compiled_cx_count": 12,
        "runtime_budget": {"external_seconds": 30, "threads": 1, "max_coordinate_points": 3, "full_family_matrix_builds": 6, "serialized_gate_list_matrix_builds": 6, "extra_prefix_and_direct_decoder_matrix_builds": 2},
        "scope": "Finite three-point matrix/serialization preflight only; no target objective, derivative, optimizer, or search.",
    }
    (HERE / "preflight_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
