"""Fit nine new run-316 seven-CNOT tail schedules with shortest linear gauges.

For each schedule, first fit arbitrary local SU(2) layers to S*T, where S is
the shortest compatible output-CNOT row gauge. Then use the fitted layers as
a start for the label-free cycle objective. This is a bounded numerical test.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup
from search13_identity_tail_fit import (
    deletion_warm, original_tail, process_fit,
)


ROOT = Path(__file__).resolve().parent
SCREEN = ROOT / "search13_all_linear_output_distance2_result.json"
OUT = ROOT / "search13_linear_gauge_nine_fit_result.json"
PRIORITY = (762, 826, 359, 740, 802, 825, 836, 847, 938)
OLD = {506, 530, 922, 943}


def cases(data):
    new = set().union(*(set(row["surviving_indices"])
                        for row in data["rows_with_survivors"])) - OLD
    assert new == set(PRIORITY), sorted(new)
    for index in PRIORITY:
        eligible = [row for row in data["rows_with_survivors"]
                    if index in row["surviving_indices"]]
        row = min(eligible, key=lambda r: (r["minimum_output_cx_length"],
                                           r["representative_output_cxs"]))
        schedule = tuple(tuple(pair) for pair in data["schedule_order"][index])
        yield index, schedule, row


def check(gates, prefix):
    candidate = {"n": 4, "gates": prefix + gates}
    result = evaluate(candidate)
    return candidate, {"cnot_count": result["cnot_count"],
                       "off_diagonal_error": result["off_diagonal_error"],
                       "valid_diagonalizer": result["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=31700)
    parser.add_argument("--process-maxiter", type=int, default=140)
    parser.add_argument("--direct-maxiter", type=int, default=200)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    data = json.loads(SCREEN.read_text())
    prefix_gates, tail_gates, old_pairs = original_tail()
    tail, prefix_matrix, _, _ = setup()
    records = []
    for index, schedule, row in cases(data):
        permutation = np.array(row["row_permutation"], dtype=int)
        assert sorted(permutation.tolist()) == list(range(16))
        target = tail[permutation, :]
        warm, deleted, distance = deletion_warm(tail_gates, old_pairs, schedule)
        rng = np.random.default_rng(args.seed + index)
        starts = [("warm", warm)]
        if index in PRIORITY[:2]:
            starts.append(("random", rng.normal(0, 0.4, (8, 4, 3))))
        trials = []
        best = None
        for kind, initial in starts:
            proc_tail, process = process_fit(schedule, initial, target,
                                             args.process_maxiter)
            _, process_check = check(proc_tail, prefix_gates)
            direct_tail, direct = direct_fit(schedule, prefix_matrix,
                                             decode_local_layers(proc_tail),
                                             args.direct_maxiter)
            candidate, direct_check = check(direct_tail, prefix_gates)
            trial = {"start": kind,
                     "process": {**process, **process_check},
                     "direct": {**direct, **direct_check}}
            trials.append(trial)
            if best is None or direct_check["off_diagonal_error"] < best[0]:
                best = direct_check["off_diagonal_error"], candidate, kind
            print(json.dumps({"index": index, "start": kind,
                              "gauge_cx_length": row["minimum_output_cx_length"],
                              "process_loss": process["loss"],
                              "direct_loss": direct["loss"],
                              "direct_off_diagonal": direct_check["off_diagonal_error"],
                              "valid": direct_check["valid_diagonalizer"]}), flush=True)
        path = ROOT / f"search13_linear_gauge_nine_fit_{index}.json"
        path.write_text(json.dumps(best[1], indent=2) + "\n")
        records.append({"index": index, "schedule": schedule,
                        "output_cxs": row["representative_output_cxs"],
                        "gauge_length": row["minimum_output_cx_length"],
                        "row_permutation": row["row_permutation"],
                        "nearest_deleted_cnot": deleted,
                        "distance_to_deletion": distance,
                        "trials": trials,
                        "best_candidate": path.name,
                        "best_start": best[2],
                        "best_off_diagonal_error": best[0]})
    result = {"source": SCREEN.name, "seed": args.seed,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "rows": records,
              "scope": "bounded numerical process and label-free fits; failure is not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"finished": len(records),
                      "best_off_diagonal_error": min(r["best_off_diagonal_error"]
                                                     for r in records)}), flush=True)


if __name__ == "__main__":
    main()
