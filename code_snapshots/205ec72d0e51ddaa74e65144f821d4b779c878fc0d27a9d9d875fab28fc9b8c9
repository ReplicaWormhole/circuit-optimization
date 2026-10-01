"""Short-fit topology tournament from run348's opening02 alternative source.

For all thirteen slots, consider the five other unordered interaction
pairs in one deterministic orientation. Exclude the unordered schedules
recorded in run347's frozen proposals/omissions, its original source, and
run340's closing-edge results. Retained proposals make two edge changes
relative to the original chain. Fit every retained proposal, then deepen
the six best fitted losses.
"""

import argparse
import hashlib
import json
from pathlib import Path

import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_fresh_chain_dynamic_rows import decode


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_alternative_opening_refine_trial0.json"
PRIOR = ROOT / "search13_alternative_opening_refine_result.json"
EXCLUDED_347 = ROOT / "search13_fresh_chain_fitted_tournament_frozen.json"
EXCLUDED_340 = ROOT / "search13_fresh_chain_closing_edge_result.json"
OUT = ROOT / "search13_opening02_fitted_tournament_result.json"
FROZEN = ROOT / "search13_opening02_fitted_tournament_frozen.json"
PAIRS = tuple((a, b) for a in range(4) for b in range(a+1, 4))
torch.set_num_threads(1)


def unordered(schedule):
    return tuple(tuple(sorted(edge)) for edge in schedule)


def checked(name, schedule, angles):
    candidate = serialize([[] for _ in range(14)], schedule, angles)
    path = ROOT / name
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short-maxiter", type=int, default=40)
    parser.add_argument("--deep-maxiter", type=int, default=400)
    parser.add_argument("--deep-count", type=int, default=6)
    args = parser.parse_args()
    assert not OUT.exists() and not FROZEN.exists()
    source = json.loads(SOURCE.read_text())
    angles = decode(source)
    source_schedule = tuple((g["control"], g["target"])
                            for g in source["gates"] if g["gate"] == "cx")
    assert len(source_schedule) == 13
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed,
                          matrix_for_schedule(source_schedule), v, False)
    expected = next(row["loss"] for row in json.loads(PRIOR.read_text())["rows"]
                    if row["trial"] == 0)
    assert abs(roundtrip - expected) < 1e-10
    old347 = json.loads(EXCLUDED_347.read_text())
    old340 = json.loads(EXCLUDED_340.read_text())
    exclusion_records = []
    for section in ("proposals", "omitted"):
        for row in old347[section]:
            exclusion_records.append({"origin_file": EXCLUDED_347.name,
                                      "section": section, "proposal": row["proposal"],
                                      "schedule": row["schedule"]})
    for row in old340["rows"]:
        exclusion_records.append({"origin_file": EXCLUDED_340.name,
                                  "slot": row["slot"], "schedule": row["schedule"]})
    exclusion_records.append({"origin_file": EXCLUDED_347.name,
                              "section": "source_schedule",
                              "schedule": old347["source_schedule"]})
    known = {unordered(row["schedule"]) for row in exclusion_records}
    assert len(known) == 66
    proposals = []
    omitted = []
    proposed = 0
    for slot, old in enumerate(source_schedule):
        for pair in PAIRS:
            if pair == tuple(sorted(old)):
                continue
            replacement = pair if slot % 2 == 0 else pair[::-1]
            schedule = list(source_schedule)
            schedule[slot] = replacement
            row = {"proposal": proposed, "slot": slot, "old": old,
                   "replacement": replacement, "schedule": schedule}
            proposed += 1
            if unordered(schedule) in known:
                omitted.append(row)
            else:
                proposals.append(row)
    assert proposed == 65 and len(omitted) == 5 and len(proposals) == 60
    assert len({unordered(row["schedule"]) for row in proposals}) == 60
    original = tuple(map(tuple, old347["source_schedule"]))
    assert all(sum(tuple(a) != tuple(b) for a, b in zip(row["schedule"], original)) == 2
               for row in proposals)
    frozen = {"source": SOURCE.name,
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "source_schedule": source_schedule,
              "excluded_files": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in (EXCLUDED_347, EXCLUDED_340)},
              "exclusion_scope": "run347 frozen proposals/omitted/original source plus run340 saved closing-edge schedules; 66 distinct unordered schedules, no broader historical claim",
              "exclusion_records": exclusion_records,
              "original_chain_schedule": original,
              "proposed": proposed, "omitted": omitted, "proposals": proposals}
    FROZEN.write_text(json.dumps(frozen, indent=2) + "\n")
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "roundtrip_error": abs(roundtrip - expected),
                      "proposed": proposed, "excluded": len(omitted),
                      "remaining": len(proposals)}), flush=True)
    short_rows = []
    fitted_angles = {}
    for proposal in proposals:
        schedule = tuple(map(tuple, proposal["schedule"]))
        fitted, record = fit(schedule, angles, fixed, v, args.short_maxiter)
        check = checked(f"search13_opening02_fitted_tournament_p{proposal['proposal']}_short.json",
                        schedule, fitted)
        row = {**proposal, **record, **check}
        short_rows.append(row)
        fitted_angles[row["proposal"]] = fitted
        print(json.dumps({"stage": "short", "proposal": row["proposal"],
                          "slot": row["slot"], "loss": row["loss"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
        if row["valid_diagonalizer"]:
            break
    deep_rows = []
    if not any(row["valid_diagonalizer"] for row in short_rows):
        selected = sorted(short_rows, key=lambda row: (row["loss"], row["proposal"]))[:args.deep_count]
        assert len({unordered(row["schedule"]) for row in selected}) == len(selected)
        for rank, parent in enumerate(selected):
            schedule = tuple(map(tuple, parent["schedule"]))
            fitted, record = fit(schedule, fitted_angles[parent["proposal"]],
                                 fixed, v, args.deep_maxiter)
            check = checked(f"search13_opening02_fitted_tournament_p{parent['proposal']}_deep.json",
                            schedule, fitted)
            row = {"rank": rank, "proposal": parent["proposal"],
                   "slot": parent["slot"], "old": parent["old"],
                   "replacement": parent["replacement"], "schedule": parent["schedule"],
                   "short_candidate": parent["candidate"],
                   "short_loss": parent["loss"], **record, **check}
            deep_rows.append(row)
            print(json.dumps({"stage": "deep", "proposal": row["proposal"],
                              "slot": row["slot"], "loss": row["loss"],
                              "off_diagonal_error": row["off_diagonal_error"],
                              "valid": row["valid_diagonalizer"]}), flush=True)
            if row["valid_diagonalizer"]:
                break
    best = min(short_rows + deep_rows, key=lambda row: row["off_diagonal_error"])
    result = {**frozen, "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "frozen_proposals": FROZEN.name,
              "initialization": "same source local angles at every short fit; deep fits start from corresponding short fitted angles; no random seeds",
              "short_maxiter": args.short_maxiter, "deep_maxiter": args.deep_maxiter,
              "short_count": len(short_rows), "deep_count": len(deep_rows),
              "short_rows": short_rows, "deep_rows": deep_rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "sixty specified one-slot mutations of opening02 source, each making two changes from original chain; one warm start, six selected by fitted loss; bounded numerical evidence only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
