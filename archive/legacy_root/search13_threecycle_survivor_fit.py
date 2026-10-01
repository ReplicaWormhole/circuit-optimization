"""Bounded fits of two three-cycle gauges on two rank-compatible tails.

Choose one gauge with the smallest modular mixed-cut rank sum among all
run-323 survivors and one with the smallest sum among three-cycles entirely
within a single eigenvalue sector. Fit pair orders 506 and 922, respectively,
first to the fixed gauged tail and then to the label-free cycle objective.
These fits are numerical search evidence, not exclusions on failure.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup
from search13_identity_tail_fit import deletion_warm, original_tail, process_fit
from search13_output_cx_near_tail_rank import MASKS, target_mod
from search13_rank8_neighborhood_mixedcut import rank_mod, realignment


ROOT = Path(__file__).resolve().parent
SCREEN = ROOT / "search13_threecycle_distance2_rank_result.json"
OUT = ROOT / "search13_threecycle_survivor_fit_result.json"
PRIME = 97


def select_gauges(data, labels):
    base = target_mod(PRIME)
    scored = []
    for row in data["rows_with_survivors"]:
        permutation = list(row["row_permutation"])
        target = base[permutation, :]
        ranks = [rank_mod(realignment(target, mask), PRIME)
                 for mask in MASKS]
        same_sector = len({labels[i] for i in row["output_cycle"]}) == 1
        scored.append((sum(ranks), sum(r == 16 for r in ranks),
                       tuple(row["output_cycle"]), same_sector, row))
    overall = min(scored, key=lambda x: x[:3])
    within = min((x for x in scored if x[3]), key=lambda x: x[:3])
    assert 506 in overall[4]["surviving_schedule_indices"]
    assert 922 in within[4]["surviving_schedule_indices"]
    return (("overall", 506, overall), ("within_sector", 922, within))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=32800)
    parser.add_argument("--process-maxiter", type=int, default=100)
    parser.add_argument("--direct-maxiter", type=int, default=180)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    data = json.loads(SCREEN.read_text())
    prefix_gates, tail_gates, original_pairs = original_tail()
    tail, prefix_matrix, labels, _ = setup()
    cases = select_gauges(data, labels)
    records = []
    for case_name, index, scored in cases:
        rank_sum, rank16, cycle, same_sector, row = scored
        schedule = tuple(tuple(edge) for edge in data["schedule_order"][index])
        target = tail[list(row["row_permutation"]), :]
        warm, deleted, distance = deletion_warm(tail_gates, original_pairs,
                                                 schedule)
        rng = np.random.default_rng(args.seed + index)
        for start_name, initial in (("deletion_warm", warm),
                                    ("random", rng.normal(0, 0.4, (8, 4, 3)))):
            process_gates, process = process_fit(schedule, initial, target,
                                                 args.process_maxiter)
            direct_gates, direct = direct_fit(
                schedule, prefix_matrix,
                decode_local_layers(process_gates), args.direct_maxiter)
            candidate = {"n": 4, "gates": prefix_gates + direct_gates}
            check = evaluate(candidate)
            path = ROOT / f"search13_threecycle_{case_name}_{start_name}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            records.append({"case": case_name, "schedule_index": index,
                            "schedule": schedule, "output_cycle": cycle,
                            "same_eigenvalue_sector": same_sector,
                            "mod97_rank_sum": rank_sum,
                            "mod97_rank16_masks": rank16,
                            "deleted_exact_tail_cnot": deleted,
                            "distance_to_deletion": distance,
                            "start": start_name, "process": process,
                            "direct": direct, "candidate": path.name,
                            "cnot_count": check["cnot_count"],
                            "off_diagonal_error": check["off_diagonal_error"],
                            "valid_diagonalizer": check["valid_diagonalizer"]})
            print(json.dumps({"case": case_name, "start": start_name,
                              "output_cycle": cycle,
                              "process_loss": process["loss"],
                              "direct_loss": direct["loss"],
                              "off_diagonal_error": check["off_diagonal_error"],
                              "valid": check["valid_diagonalizer"]}), flush=True)
    best = min(records, key=lambda r: r["off_diagonal_error"])
    result = {"screen": SCREEN.name, "seed": args.seed,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "rows": records,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "two rank-profile-selected three-cycle gauges, two seven-tail orders and two starts each; numerical failure is not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
