"""Numerical nullspace structure of the refined simultaneous rank-eight gauge.

Inspect left/right rank-eight kernel projectors at mixed masks 83 and 163,
and amplitude patterns of the resulting output gauge. This is a diagnostic
for an exact sparse ansatz, not a rank certificate.
"""

import json
from pathlib import Path

import numpy as np
import torch

from search13_fulltail_invariant_gauge_rank import gauge, setup
from search13_rank8_bridge import exact_gauge


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_mixedcut_gauge_refine_result.json"
OUT = ROOT / "search13_joint_mixedcut_kernel_ansatz_result.json"
MASKS = (83, 163)


def realign(matrix, mask):
    left = tuple(j for j in range(8) if mask & (1 << j))
    right = tuple(j for j in range(8) if j not in left)
    return matrix.reshape((2,) * 8).transpose(left + right).reshape(16, 16)


def entries(matrix, threshold):
    return [[int(i), int(j), float(z.real), float(z.imag)]
            for i, row in enumerate(matrix) for j, z in enumerate(row)
            if abs(z) > threshold]


def main():
    assert not OUT.exists()
    data = json.loads(SOURCE.read_text())
    tail, _, _, sectors = setup()
    with torch.no_grad():
        h = gauge(torch.tensor(data["refined_angles"], dtype=torch.float64),
                  sectors).numpy()
    composed = h @ exact_gauge()
    target = composed @ tail
    cuts = {}
    for mask in MASKS:
        matrix = realign(target, mask)
        u, singular, vh = np.linalg.svd(matrix)
        right = vh[8:].conj().T @ vh[8:]
        left = u[:, 8:] @ u[:, 8:].conj().T
        cuts[str(mask)] = {
            "singular_values": singular.tolist(),
            "right_kernel_projector_diagonal": np.real(np.diag(right)).tolist(),
            "left_kernel_projector_diagonal": np.real(np.diag(left)).tolist(),
            "right_kernel_projector_entries_above_0_2": entries(right, 0.2),
            "left_kernel_projector_entries_above_0_2": entries(left, 0.2),
            "near_zero_columns": [j for j in range(16)
                                  if np.linalg.norm(matrix[:, j]) < 1e-8],
            "near_zero_rows": [i for i in range(16)
                               if np.linalg.norm(matrix[i]) < 1e-8],
        }
    row_magnitudes = [[float(abs(z)) for z in row] for row in composed]
    result = {"source": SOURCE.name,
              "output_gauge_row_magnitudes": row_magnitudes,
              "output_gauge_entries_above_0_1": entries(composed, 0.1),
              "cuts": cuts,
              "scope": "numerical SVD/projector diagnostics only; singular values and near-zero relations require exactification"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"mask_singular_ninth": {str(k): cuts[str(k)]["singular_values"][8]
                                             for k in MASKS},
                      "gauge_entries_above_0_1": len(result["output_gauge_entries_above_0_1"]),
                      "near_zero_columns": {str(k): cuts[str(k)]["near_zero_columns"] for k in MASKS}}))


if __name__ == "__main__":
    main()
