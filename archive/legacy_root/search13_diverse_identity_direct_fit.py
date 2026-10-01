"""Diverse identity-tail schedules: direct cycle fits with two-stage budget.

Select three schedules from each balanced crossing profile, favoring large
pairwise Hamming distance and a spread of edit distances from the exact tail.
For every selected schedule, fit the label-free cycle objective from deletion
warm and random local SU(2) starts. Deepen six strongest distinct-profile
trials. This is a numerical search, not a lower-bound proof.
"""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup
from search13_identity_tail_fit import deletion_warm, original_tail


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_identity_tail_topology_result.json"
OUT = ROOT / "search13_diverse_identity_direct_fit_result.json"


def schedule(row):
    return tuple(tuple(pair) for pair in row["schedule"])


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def select(rows, per_profile=3):
    eligible = [r for r in rows if r["distance_to_one_exact8_deletion"] >= 3]
    profiles = sorted({tuple(r["balanced_crossings"]) for r in eligible})
    assert len(profiles) == 6
    selected = []
    for slot in range(per_profile):
        for profile in profiles:
            pool = [r for r in eligible
                    if tuple(r["balanced_crossings"]) == profile and r not in selected]
            assert pool
            # Spread edit distances as well as ordered pair schedules.
            target_distance = (3, 5, 7)[slot]
            def score(row):
                s = schedule(row)
                spacing = min((hamming(s, schedule(old)) for old in selected),
                              default=7)
                same_profile_spacing = min((hamming(s, schedule(old))
                                            for old in selected
                                            if tuple(old["balanced_crossings"]) == profile),
                                           default=7)
                return (same_profile_spacing, spacing,
                        -abs(row["distance_to_one_exact8_deletion"] - target_distance),
                        -row["adjacent_repeats"], s)
            selected.append(max(pool, key=score))
    assert len(selected) == 18
    assert len({schedule(r) for r in selected}) == 18
    return selected


def checked(prefix, tail):
    candidate = {"n": 4, "gates": prefix + tail}
    result = evaluate(candidate)
    return candidate, {"cnot_count": result["cnot_count"],
                       "off_diagonal_error": result["off_diagonal_error"],
                       "valid_diagonalizer": result["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=32100)
    parser.add_argument("--screen-maxiter", type=int, default=70)
    parser.add_argument("--deepen-maxiter", type=int, default=350)
    parser.add_argument("--deepen-count", type=int, default=6)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    assert len(source["surviving_schedules"]) == 1428
    rows = select(source["surviving_schedules"])
    prefix_gates, tail_gates, original_pairs = original_tail()
    _, prefix_matrix, _, _ = setup()
    records = []
    best_candidates = {}
    for case, row in enumerate(rows):
        edges = schedule(row)
        warm, deleted, distance = deletion_warm(tail_gates, original_pairs, edges)
        assert distance == row["distance_to_one_exact8_deletion"]
        rng = np.random.default_rng(args.seed + 137*case)
        starts = (("structured_deletion", warm),
                  ("random", rng.normal(0, 0.4, (8, 4, 3))))
        trials = []
        for kind, initial in starts:
            fitted_tail, fit = direct_fit(edges, prefix_matrix, initial,
                                          args.screen_maxiter)
            candidate, check = checked(prefix_gates, fitted_tail)
            trials.append({"start": kind, "fit": fit, "checker": check})
            if case not in best_candidates or fit["loss"] < best_candidates[case][0]:
                best_candidates[case] = (fit["loss"], candidate, fitted_tail, kind)
            print(json.dumps({"case": case, "profile": row["balanced_crossings"],
                              "distance": distance, "start": kind,
                              "loss": fit["loss"],
                              "off_diagonal": check["off_diagonal_error"],
                              "valid": check["valid_diagonalizer"]}), flush=True)
        records.append({"case": case, "schedule": edges,
                        "balanced_crossings": row["balanced_crossings"],
                        "distance_to_one_exact8_deletion": distance,
                        "adjacent_repeats": row["adjacent_repeats"],
                        "nearest_deleted_cnot": deleted,
                        "screen_trials": trials})

    # Deepen the strongest four plus distinct-profile alternatives when useful.
    ranked = sorted(range(len(records)), key=lambda case: best_candidates[case][0])
    deep = ranked[:min(4, args.deepen_count)]
    represented = {tuple(records[case]["balanced_crossings"]) for case in deep}
    for case in ranked:
        profile = tuple(records[case]["balanced_crossings"])
        if len(deep) < args.deepen_count and profile not in represented:
            deep.append(case)
            represented.add(profile)
    for case in ranked:
        if len(deep) < args.deepen_count and case not in deep:
            deep.append(case)
    for case in deep:
        edges = schedule(rows[case])
        _, _, fitted_tail, start_kind = best_candidates[case]
        deeper_tail, fit = direct_fit(edges, prefix_matrix,
                                      decode_local_layers(fitted_tail),
                                      args.deepen_maxiter)
        candidate, check = checked(prefix_gates, deeper_tail)
        records[case]["deepen"] = {"from_start": start_kind,
                                   "fit": fit, "checker": check}
        if fit["loss"] < best_candidates[case][0]:
            best_candidates[case] = (fit["loss"], candidate, deeper_tail,
                                     "deepened_" + start_kind)
        print(json.dumps({"deepen_case": case,
                          "loss": fit["loss"],
                          "off_diagonal": check["off_diagonal_error"],
                          "valid": check["valid_diagonalizer"]}), flush=True)

    for case, record in enumerate(records):
        loss, candidate, _, kind = best_candidates[case]
        path = ROOT / f"search13_diverse_identity_direct_fit_{case:02d}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        final = evaluate(candidate)
        record["best_candidate"] = path.name
        record["best_start"] = kind
        record["best_loss"] = loss
        record["best_off_diagonal_error"] = final["off_diagonal_error"]
        record["best_valid_diagonalizer"] = final["valid_diagonalizer"]
    best = min(records, key=lambda r: r["best_off_diagonal_error"])
    result = {"source": SOURCE.name, "seed": args.seed,
              "screen_maxiter": args.screen_maxiter,
              "deepen_maxiter": args.deepen_maxiter,
              "deepen_cases": deep,
              "selection": "three diverse schedules from each of six balanced crossing profiles, edit distance >=3",
              "selected_distance_counts": dict(Counter(r["distance_to_one_exact8_deletion"]
                                               for r in records)),
              "rows": records,
              "best_checker": {"case": best["case"],
                               "candidate": best["best_candidate"],
                               "off_diagonal_error": best["best_off_diagonal_error"],
                               "valid_diagonalizer": best["best_valid_diagonalizer"]},
              "scope": "bounded gauge-free direct numerical fits; failures are not no-go results"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"finished": len(records),
                      "deepen": len(deep),
                      "best": result["best_checker"]}), flush=True)


if __name__ == "__main__":
    main()
