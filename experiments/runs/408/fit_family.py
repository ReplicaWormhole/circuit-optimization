"""Full local SU2 freedom around an arbitrary eleven-CX chronological word.

Importing this module constructs no circuit matrices, tensors or random draws.
"""
from __future__ import annotations
import math
import numpy as np
import torch


def cx_numpy(control, target):
    out = np.zeros((16, 16), dtype=np.complex128)
    for col in range(16):
        row = col ^ (1 << (3-target)) if col & (1 << (3-control)) else col
        out[row, col] = 1
    return out


def right_shift_numpy():
    out = np.zeros((16, 16), dtype=np.complex128)
    for col in range(16):
        out[((col & 1) << 3) | (col >> 1), col] = 1
    return out


def rotvec_numpy(r):
    r = np.asarray(r, dtype=float)
    length = float(np.linalg.norm(r))
    factor = math.sin(length/2)/length if length else 0.5
    x, y, z = r * factor
    c = math.cos(length/2)
    return np.array([[c-1j*z, -y-1j*x], [y-1j*x, c+1j*z]], dtype=np.complex128)


def matrix_rotvec(U):
    """SU2 principal rotation vector; discarded determinant phase is global."""
    U = np.asarray(U, dtype=np.complex128)
    Q = U / np.sqrt(np.linalg.det(U))
    if float(np.trace(Q).real) < 0:
        Q = -Q
    scalar = float(np.trace(Q).real / 2)
    vector = np.array([-((Q[0, 1]+Q[1, 0]).imag)/2,
                       ((Q[1, 0]-Q[0, 1]).real)/2,
                       -((Q[0, 0]-Q[1, 1]).imag)/2])
    length = float(np.linalg.norm(vector))
    if length < 1e-14:
        return np.zeros(3)
    return vector * (2 * math.atan2(length, scalar) / length)


def collapse(candidate, delete_cx_ordinal, checker):
    """Delete one CX and merge each wire's local gates in each free layer."""
    layers = []
    current = [np.eye(2, dtype=np.complex128) for q in range(4)]
    topology = []
    deleted = []
    cx_ordinal = 0
    deleted_source_position = None
    for position, gate in enumerate(candidate['gates']):
        if gate['gate'] == 'cx':
            if cx_ordinal == delete_cx_ordinal:
                deleted_source_position = position
            else:
                layers.append(current)
                current = [np.eye(2, dtype=np.complex128) for q in range(4)]
                topology.append([gate['control'], gate['target']])
                deleted.append(dict(gate))
            cx_ordinal += 1
        else:
            q = gate['qubit']
            current[q] = checker.one_qubit_matrix(gate) @ current[q]
            deleted.append(dict(gate))
    layers.append(current)
    if cx_ordinal != 12 or deleted_source_position is None or len(topology) != 11 or len(layers) != 12:
        raise ValueError('expected exactly one of twelve source CX to be deleted')
    vectors = np.array([[matrix_rotvec(U) for U in layer] for layer in layers]).reshape(144)
    return vectors, topology, {'n': 4, 'gates': deleted}, layers, deleted_source_position


def unitary_numpy(vectors, topology):
    layers = np.asarray(vectors).reshape(12, 4, 3)
    U = np.eye(16, dtype=np.complex128)
    for index, layer in enumerate(layers):
        M = rotvec_numpy(layer[0])
        for q in range(1, 4):
            M = np.kron(M, rotvec_numpy(layer[q]))
        U = M @ U
        if index < len(topology):
            U = cx_numpy(*topology[index]) @ U
    return U


def unitary_from_gate_list(candidate, checker):
    U = np.eye(16, dtype=np.complex128)
    for gate in candidate['gates']:
        if gate['gate'] == 'cx':
            M = cx_numpy(gate['control'], gate['target'])
        else:
            M = checker.embedded_one_qubit(checker.one_qubit_matrix(gate), gate['qubit'], 4)
        U = M @ U
    return U


def global_phase_error(U, V):
    overlap = np.vdot(V, U)
    phase = overlap/abs(overlap) if abs(overlap) else 1
    return float(np.max(np.abs(U-phase*V)))


def u3_from_su2(U, qubit):
    a, b, c, d = U[0, 0], U[0, 1], U[1, 0], U[1, 1]
    theta = 2*math.atan2(abs(c), abs(a))
    if abs(a) > 1e-12:
        phase = np.angle(a)
        phi = np.angle(c)-phase if abs(c) > 1e-12 else 0.0
        lam = np.angle(-b)-phase if abs(b) > 1e-12 else np.angle(d)-phase
    else:
        phase = np.angle(c)
        phi = 0.0
        lam = np.angle(-b)-phase
    return {'gate': 'u3', 'qubit': qubit, 'theta': float(theta), 'phi': float(phi), 'lam': float(lam)}


def serialize(vectors, topology):
    layers = np.asarray(vectors).reshape(12, 4, 3)
    gates = []
    for index, layer in enumerate(layers):
        gates.extend(u3_from_su2(rotvec_numpy(layer[q]), q) for q in range(4))
        if index < len(topology):
            gates.append({'gate': 'cx', 'control': topology[index][0], 'target': topology[index][1]})
    return {'n': 4, 'gates': gates, 'metadata': {
        'description': 'eleven fixed CX with all144 SU2 rotation-vector coordinates free',
        'evidence_scope': 'bounded numerical refit; not an exact certificate'}}


class TorchWord:
    def __init__(self, topology):
        self.pauli = torch.tensor([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], dtype=torch.complex128)
        self.eye2 = torch.eye(2, dtype=torch.complex128)
        self.cx = [torch.tensor(cx_numpy(*edge), dtype=torch.complex128) for edge in topology]
        self.target = torch.tensor(right_shift_numpy(), dtype=torch.complex128)

    def unitary(self, x):
        r = x.reshape(12, 4, 3)
        lengths = torch.linalg.vector_norm(r, dim=-1)
        # torch.sinc(t)=sin(pi*t)/(pi*t); this limit is regular at zero.
        factors = 0.5 * torch.sinc(lengths/(2*math.pi))
        axes = torch.sum(r[..., None, None] * self.pauli, dim=-3)
        locals_ = torch.cos(lengths/2)[..., None, None]*self.eye2 - 1j*factors[..., None, None]*axes
        U = torch.eye(16, dtype=torch.complex128)
        for index, layer in enumerate(locals_):
            M = layer[0].contiguous()
            for q in range(1, 4):
                M = torch.kron(M, layer[q].contiguous())
            U = M @ U
            if index < len(self.cx):
                U = self.cx[index] @ U
        return U

    def loss(self, x):
        U = self.unitary(x)
        A = U @ self.target @ U.conj().T
        off = A - torch.diag(torch.diagonal(A))
        return torch.sum(torch.abs(off)**2).real/16
