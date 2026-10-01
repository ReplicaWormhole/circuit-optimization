"""Fit a constrained exact-angle four-matchgate tail with local Z phases.

First disjoint pair layer (02,13): Weyl (pi/4,pi/8,0).
Second disjoint pair layer (01,23): Weyl (pi/8,pi/8,0).
Use only Z rotations before and between layers after the exact six-CX
prefix; output Z rotations cannot affect V4 diagonalization. Search four
relative YY sign patterns and two random starts with 500 Adam steps each.
This is a restricted eight-parameter model, not general exactification.
"""

import json
import math
from pathlib import Path

import numpy as np
import torch

from check_circuit import (cnot, embedded_one_qubit, one_qubit_matrix,
                           right_shift)


ROOT = Path(__file__).parent
torch.set_num_threads(1)
DTYPE = torch.complex128
I = torch.eye(16, dtype=DTYPE)
Z2 = torch.diag(torch.tensor([1.0, -1.0], dtype=DTYPE))
X2 = torch.tensor([[0, 1], [1, 0]], dtype=DTYPE)
Y2 = torch.tensor([[0, -1j], [1j, 0]], dtype=DTYPE)


def embed_two(pauli, q, r):
    result = torch.ones((1, 1), dtype=DTYPE)
    for j in range(4):
        result = torch.kron(result, pauli if j in (q, r)
                            else torch.eye(2, dtype=DTYPE))
    return result


def entangler(pair, a, b):
    xx = embed_two(X2, *pair)
    yy = embed_two(Y2, *pair)
    return ((math.cos(a) * I - 1j * math.sin(a) * xx) @
            (math.cos(b) * I - 1j * math.sin(b) * yy))


def local_z(angles):
    result = torch.ones((1, 1), dtype=DTYPE)
    for angle in angles:
        local = torch.cos(angle / 2) * torch.eye(2, dtype=DTYPE) - (
            1j * torch.sin(angle / 2) * Z2)
        result = torch.kron(result, local)
    return result


def prefix_matrix():
    candidate = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    matrix = np.eye(16, dtype=complex)
    for gate in candidate["gates"][:13]:
        fixed = (cnot(gate["control"], gate["target"], 4)
                 if gate["gate"] == "cx"
                 else embedded_one_qubit(
                     one_qubit_matrix(gate), gate["qubit"], 4))
        matrix = fixed @ matrix
    return torch.tensor(matrix, dtype=DTYPE)


def fit(sign_first, sign_second, seed, prefix, v):
    first = (entangler((0, 2), math.pi / 4,
                       sign_first * math.pi / 8) @
             entangler((1, 3), math.pi / 4,
                       sign_first * math.pi / 8))
    second = (entangler((0, 1), math.pi / 8,
                        sign_second * math.pi / 8) @
              entangler((2, 3), math.pi / 8,
                        sign_second * math.pi / 8))
    generator = torch.Generator().manual_seed(seed)
    angles = torch.nn.Parameter(
        0.2 * torch.randn(2, 4, generator=generator,
                          dtype=torch.float64))
    optimizer = torch.optim.Adam([angles], lr=0.04)
    best = float("inf")
    best_angles = None
    for _ in range(500):
        optimizer.zero_grad()
        u = second @ local_z(angles[1]) @ first @ (
            local_z(angles[0]) @ prefix)
        conjugated = u @ v @ u.conj().T
        offdiag = conjugated - torch.diag(torch.diag(conjugated))
        objective = torch.mean(torch.abs(offdiag) ** 2)
        value = float(objective.detach())
        if value < best:
            best = value
            best_angles = angles.detach().clone().numpy()
        objective.backward()
        optimizer.step()
    return {"sign_first": sign_first, "sign_second": sign_second,
            "seed": seed, "best_loss": best,
            "angles_over_pi": (best_angles / math.pi).tolist()}


def main():
    prefix = prefix_matrix()
    v = torch.tensor(right_shift(4), dtype=DTYPE)
    records = [fit(sf, ss, seed, prefix, v)
               for sf in (-1, 1) for ss in (-1, 1)
               for seed in (0, 1)]
    records.sort(key=lambda row: row["best_loss"])
    result = {"trials": len(records), "steps_each": 500,
              "best": records[0], "records": records}
    (ROOT / "algebraic14_matchgate_z_model_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({"trials": len(records), "best": records[0]},
                     indent=2))


if __name__ == "__main__":
    main()
