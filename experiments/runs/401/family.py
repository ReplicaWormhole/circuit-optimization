"""Differentiable ordered full98 family for the Bell-decoder word."""
import importlib.util
import math
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from check_circuit import cnot, right_shift, evaluate as circuit_evaluate  # noqa: E402
from native_gate_check import evaluate as native_evaluate  # noqa: E402
from delete14_search import axis_rotation, axis_to_u3, kron_gate  # noqa: E402

HELPER_PATH = ROOT / "collaboration/work/verifier_label12/parity12/native_helpers.py"
_spec = importlib.util.spec_from_file_location("bellprefix_native_helpers", HELPER_PATH)
helper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(helper)

SCHEDULE = [
    ("xx_yy", 0, 2),
    ("xx_yy", 1, 3),
    ("cx", 0, 1),
    ("cx", 2, 1),
    ("cx", 3, 1),
    ("cx", 1, 2),
    ("xx_yy", 0, 1),
    ("xx_yy", 2, 3),
]
F_PAIRS = [(0, 2), (1, 3), (0, 1), (2, 3)]
CX_PAIRS = [(0, 1), (2, 1), (3, 1), (1, 2)]


def matrix(z):
    """Matrix of L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6."""
    local = z[:84].reshape(7, 4, 3)
    ab = z[84:90].reshape(2, 3)
    f = z[90:].reshape(4, 2)
    u = torch.eye(16, dtype=torch.complex128)

    def apply_layer(layer):
        nonlocal u
        for q in range(4):
            u = kron_gate(axis_rotation(*local[layer, q]), q) @ u

    def apply_cx(index):
        nonlocal u
        control, target = CX_PAIRS[index]
        u = torch.as_tensor(cnot(control, target, 4), dtype=torch.complex128) @ u

    apply_layer(0)
    u = helper.f_matrix(f[0, 0], f[0, 1], *F_PAIRS[0]) @ u
    apply_layer(1)
    u = helper.f_matrix(f[1, 0], f[1, 1], *F_PAIRS[1]) @ u
    apply_layer(2)
    apply_cx(0)
    u = kron_gate(axis_rotation(*ab[0]), 1) @ u
    apply_cx(1)
    u = kron_gate(axis_rotation(*ab[1]), 1) @ u
    apply_cx(2)
    apply_layer(3)
    apply_cx(3)
    apply_layer(4)
    u = helper.f_matrix(f[2, 0], f[2, 1], *F_PAIRS[2]) @ u
    apply_layer(5)
    u = helper.f_matrix(f[3, 0], f[3, 1], *F_PAIRS[3]) @ u
    apply_layer(6)
    return u


def objective(z):
    u = matrix(z)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    b = u @ v @ u.conj().T
    off = b - torch.diag(torch.diagonal(b))
    return (off.abs() ** 2).sum().real / 16


def serialize(z):
    local = np.asarray(z[:84]).reshape(7, 4, 3)
    ab = np.asarray(z[84:90]).reshape(2, 3)
    f = np.asarray(z[90:]).reshape(4, 2)
    gates = []

    def add_layer(k):
        for q in range(4):
            gates.append({"gate": "u3", "qubit": q, **axis_to_u3(*local[k, q])})

    def add_cx(index):
        control, target = CX_PAIRS[index]
        gates.append({"gate": "cx", "control": control, "target": target})

    def add_f(index):
        gates.append(helper.native_gate(*F_PAIRS[index], *f[index]))

    add_layer(0)
    add_f(0)
    add_layer(1)
    add_f(1)
    add_layer(2)
    add_cx(0)
    gates.append({"gate": "u3", "qubit": 1, **axis_to_u3(*ab[0])})
    add_cx(1)
    gates.append({"gate": "u3", "qubit": 1, **axis_to_u3(*ab[1])})
    add_cx(2)
    add_layer(3)
    add_cx(3)
    add_layer(4)
    add_f(2)
    add_layer(5)
    add_f(3)
    add_layer(6)
    return {"n": 4, "gates": gates}


def bell_decoder_point():
    """Exact SO(3) chart point for E02^dagger E13^dagger prefix."""
    x = np.zeros(98, dtype=np.float64)
    local = x[:84].reshape(7, 4, 3)
    f = x[90:].reshape(4, 2)
    # Ry(pi/2) S^dagger is, up to phase, the SO(3) rotation vector below.
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


def initial_point(seed, sigma):
    rng = np.random.default_rng(seed)
    base = bell_decoder_point()
    noise = rng.normal(0.0, sigma, size=98)
    return base, noise, base + noise


def compile_native(candidate):
    return helper.compile_native(candidate)


def evaluate_native(candidate):
    return native_evaluate(candidate)


def evaluate_compiled(candidate):
    return circuit_evaluate(candidate)
