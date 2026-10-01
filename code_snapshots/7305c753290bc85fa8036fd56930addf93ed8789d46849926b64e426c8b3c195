"""One-step directed-CNOT topology escape from the refined run-327 basin.

Score every replacement of exactly one of the thirteen directed CNOT edges
while retaining the 168 optimized local angles. Fit the six strongest raw
mutations on distinct slots that change the unordered pair. All local layers
are free during each fit. This is bounded numerical evidence only.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import (
    fit, matrix_for_schedule, objective, serialize,
)
from search13_joint_prefix_orbit_fit import CASES
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_joint_prefix_orbit_refine_0.json"
OUT = ROOT / "search13_joint_prefix_one_step_escape_result.json"
PAIRS = tuple((a, b) for a in range(4) for b in range(4) if a != b)
torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=int, default=6)
    parser.add_argument("--maxiter", type=int, default=250)
    args = parser.parse_args()
    source = json.loads(SOURCE.read_text())
    angles = decode_layers(source)
    deleted, rewire, blocks = CASES[1]
    chunks, old_cxs, fixed, _ = prepare(deleted, blocks,
                                           prefix_rewire=rewire)
    schedule = tuple((gate["control"], gate["target"])
                     for gate in old_cxs)
    assert sum(g["gate"] == "cx" for g in source["gates"]) == 13
    assert tuple((g["control"], g["target"])
                 for g in source["gates"] if g["gate"] == "cx") == schedule
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    baseline = objective(angles.ravel(), fixed,
                         matrix_for_schedule(schedule), v, False)
    assert abs(baseline - 0.0787330125161619) < 1e-8, baseline

    screen = []
    for slot, old in enumerate(schedule):
        for replacement in PAIRS:
            if replacement == old:
                continue
            mutated = list(schedule)
            mutated[slot] = replacement
            raw = objective(angles.ravel(), fixed,
                            matrix_for_schedule(mutated), v, False)
            screen.append({"slot": slot, "old": old, "replacement": replacement,
                           "same_unordered_pair": set(replacement) == set(old),
                           "raw_loss": raw})
    assert len(screen) == 13 * 11
    eligible = sorted((row for row in screen if not row["same_unordered_pair"]),
                      key=lambda row: (row["raw_loss"], row["slot"],
                                       row["replacement"]))
    selected = []
    used_slots = set()
    for row in eligible:
        if row["slot"] not in used_slots:
            selected.append(row)
            used_slots.add(row["slot"])
            if len(selected) == args.fits:
                break

    rows = []
    for rank, proposal in enumerate(selected):
        mutated = list(schedule)
        mutated[proposal["slot"]] = tuple(proposal["replacement"])
        fitted, record = fit(mutated, angles, fixed, v, args.maxiter)
        candidate = serialize(chunks, mutated, fitted)
        path = ROOT / f"search13_joint_prefix_one_step_escape_rank{rank}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"rank": rank, **proposal, **record,
               "candidate": path.name, "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps({"rank": rank, "slot": row["slot"],
                          "replacement": row["replacement"],
                          "raw_loss": row["raw_loss"],
                          "loss": row["loss"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["off_diagonal_error"])
    result = {"source": SOURCE.name, "baseline_loss": baseline,
              "source_topology": schedule,
              "screened_directed_replacements": len(screen),
              "screen": screen, "selection_rule": "lowest raw loss at distinct slots, excluding pure CNOT direction reversals",
              "fits": args.fits, "maxiter": args.maxiter, "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "one-step mutations of one fixed 13-CNOT topology and warm local layers; bounded fits, no no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
