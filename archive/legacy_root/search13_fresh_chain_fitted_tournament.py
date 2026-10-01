"""Short-fit topology tournament from the fresh chain13 numerical source.

For all thirteen slots, consider the five other unordered interaction
pairs in one deterministic orientation. Exclude only the thirteen known
closing-edge schedules of run340, comparing unordered schedules. Fit every
remaining proposal briefly, then deepen the six best fitted losses.
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
SOURCE = ROOT / "search13_fresh_chain_refine_trial1.json"
PRIOR = ROOT / "search13_fresh_chain_refine_result.json"
EXCLUDED = ROOT / "search13_fresh_chain_closing_edge_result.json"
OUT = ROOT / "search13_fresh_chain_fitted_tournament_result.json"
FROZEN = ROOT / "search13_fresh_chain_fitted_tournament_frozen.json"
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
                    if row["trial"] == 1)
    assert abs(roundtrip - expected) < 1e-10
    known = {unordered(row["schedule"])
             for row in json.loads(EXCLUDED.read_text())["rows"]}
    assert len(known) == 13
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
    assert proposed == 65 and len(omitted) == 13 and len(proposals) == 52
    assert len({unordered(row["schedule"]) for row in proposals}) == 52
    frozen = {"source": SOURCE.name,
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "source_schedule": source_schedule,
              "excluded_file": EXCLUDED.name,
              "excluded_sha256": hashlib.sha256(EXCLUDED.read_bytes()).hexdigest(),
              "exclusion_scope": "only thirteen saved run340 closing-edge schedules, compared as unordered interaction schedules",
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
        check = checked(f"search13_fresh_chain_fitted_tournament_p{proposal['proposal']}_short.json",
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
            check = checked(f"search13_fresh_chain_fitted_tournament_p{parent['proposal']}_deep.json",
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
              "scope": "specified one-slot schedules, one source warm start per short fit; six selected by fitted loss; bounded numerical evidence only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
