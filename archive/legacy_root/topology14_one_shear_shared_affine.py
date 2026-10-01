"""Finite GL(3,2) search for a shared normal form of one-shear sectors.

Conjugate both exact affine 3-bit V4 sector actions by every linear map
reachable via 3-bit CNOTs (all 168 GL3 maps). Record shortest CNOT length
of the shared conjugator, unconditional z=0 linear action, and conditional
linear discrepancy T1*T0^-1. X translations are tracked separately.
"""

import argparse
import collections
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IDENTITY = tuple(range(8))
CX_EDGES = tuple((c, t) for c in range(3) for t in range(3) if c != t)


def compose(left, right):
    return tuple(left[right[x]] for x in range(8))


def inverse(permutation):
    return tuple(permutation.index(x) for x in range(8))


def cnot(c, t):
    c_mask, t_mask = 1 << (2 - c), 1 << (2 - t)
    return tuple(x ^ t_mask if x & c_mask else x for x in range(8))


def linear_bfs():
    paths = {IDENTITY: []}
    queue = collections.deque([IDENTITY])
    while queue:
        current = queue.popleft()
        for edge in CX_EDGES:
            next_perm = compose(cnot(*edge), current)
            if next_perm not in paths:
                paths[next_perm] = paths[current] + [edge]
                queue.append(next_perm)
    assert len(paths) == 168
    return paths


def affine_parts(permutation):
    b = permutation[0]
    linear = tuple(permutation[x] ^ b for x in range(8))
    assert all(linear[x ^ y] == linear[x] ^ linear[y]
               for x in range(8) for y in range(8))
    return linear, b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    action = json.loads((ROOT / "topology14_one_shear_action_result.json").read_text())
    t0 = tuple(action["sectors"]["0"]["permutation"])
    t1 = tuple(action["sectors"]["1"]["permutation"])
    paths = linear_bfs()
    rows = []
    for l, path in paths.items():
        li = inverse(l)
        a0 = compose(l, compose(t0, li))
        a1 = compose(l, compose(t1, li))
        correction = compose(a1, inverse(a0))
        a0linear, a0shift = affine_parts(a0)
        a1linear, a1shift = affine_parts(a1)
        clinear, cshift = affine_parts(correction)
        rows.append({"conjugator": path, "conjugator_cx": len(path),
                     "z0_linear_cx": len(paths[a0linear]),
                     "z0_linear_path": paths[a0linear], "z0_shift": a0shift,
                     "z1_linear_cx": len(paths[a1linear]),
                     "z1_shift": a1shift,
                     "controlled_correction_linear_cx": len(paths[clinear]),
                     "controlled_correction_linear_path": paths[clinear],
                     "controlled_correction_shift": cshift,
                     "z0_permutation": a0, "z1_permutation": a1,
                     "controlled_correction_permutation": correction})
    rows.sort(key=lambda row: (row["controlled_correction_linear_cx"],
                               row["conjugator_cx"], row["z0_linear_cx"],
                               row["controlled_correction_shift"].bit_count()))
    histogram = {}
    for row in rows:
        key = str(row["controlled_correction_linear_cx"])
        histogram[key] = histogram.get(key, 0) + 1
    result = {"gl3_count": len(paths),
              "controlled_correction_linear_cx_histogram": histogram,
              "best": rows[:12],
              "search_scope": "linear shared conjugators only; X translations tracked but not costed",
              "interpretation": "conditional discrepancy still needs parity-controlled gates; counts are affine CNOT word lengths, not full quantum circuit CX counts"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"gl3_count": len(paths),
                      "histogram": histogram,
                      "best": rows[:3]}, indent=2))


if __name__ == "__main__":
    main()
