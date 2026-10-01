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
OUT = ROOT / "search13_opening02_recovery_result.json"
STATE = ROOT / "search13_opening02_recovery_checkpoint.json"
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
    assert not OUT.exists()
    frozen = json.loads(FROZEN.read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == frozen["source_sha256"]
    source = json.loads(SOURCE.read_text())
    angles = decode(source)
    source_schedule = tuple(map(tuple, frozen["source_schedule"]))
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed, matrix_for_schedule(source_schedule), v, False)
    expected = next(row["loss"] for row in json.loads(PRIOR.read_text())["rows"] if row["trial"] == 0)
    assert abs(roundtrip - expected) < 1e-10
    proposals = frozen["proposals"]
    assert len(proposals) == 60
    recovered = json.loads((ROOT / "search13_opening02_fitted_tournament_interrupted_validation.json").read_text())["reconstructed_rows"]
    recovered_map = {row["proposal"]:row for row in recovered}
    assert len(recovered_map) == 7
    def checkpoint(short_rows, deep_rows):
        temporary = STATE.with_suffix(".tmp")
        temporary.write_text(json.dumps({"short_rows":short_rows,"deep_rows":deep_rows},indent=2)+"\n")
        temporary.replace(STATE)
    saved = json.loads(STATE.read_text()) if STATE.exists() else {"short_rows":[],"deep_rows":[]}
    saved_short = {row["proposal"]:row for row in saved["short_rows"]}
    print(json.dumps({"source_roundtrip_loss":roundtrip,"roundtrip_error":abs(roundtrip-expected),"recovered":7,"remaining_new_short":53}),flush=True)
    short_rows = []
    fitted_angles = {}
    for proposal in proposals:
        schedule = tuple(map(tuple, proposal["schedule"]))
        pid = proposal["proposal"]
        if pid in saved_short:
            row = saved_short[pid]
            fitted = decode(json.loads((ROOT / row["candidate"]).read_text()))
        elif pid in recovered_map:
            old = recovered_map[pid]
            path = ROOT / old["candidate"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == old["sha256"]
            fitted = decode(json.loads(path.read_text()))
            loss = objective(fitted.ravel(),fixed,matrix_for_schedule(schedule),v,False)
            assert abs(loss-old["loss"])<1e-11
            check = old["checker"]
            row = {**proposal,"candidate":path.name,"loss":loss,"cnot_count":check["cnot_count"],"off_diagonal_error":check["off_diagonal_error"],"valid_diagonalizer":check["valid_diagonalizer"],"recovered_from_run":351,"candidate_sha256":old["sha256"],"optimizer_metadata":"unavailable; loss reconstructed from saved gate list; not a new fit"}
        else:
            fitted, record = fit(schedule, angles, fixed, v, args.short_maxiter)
            check = checked(f"search13_opening02_recovery_p{pid}_short.json",schedule,fitted)
            row = {**proposal,**record,**check,"new_fit":True}
        short_rows.append(row)
        fitted_angles[pid] = fitted
        checkpoint(short_rows, saved["deep_rows"])
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
            if rank < len(saved["deep_rows"]):
                row = saved["deep_rows"][rank]
                assert row["proposal"] == parent["proposal"]
                deep_rows.append(row)
                continue
            fitted, record = fit(schedule, fitted_angles[parent["proposal"]],
                                 fixed, v, args.deep_maxiter)
            check = checked(f"search13_opening02_recovery_p{parent['proposal']}_deep.json",
                            schedule, fitted)
            row = {"rank": rank, "proposal": parent["proposal"],
                   "slot": parent["slot"], "old": parent["old"],
                   "replacement": parent["replacement"], "schedule": parent["schedule"],
                   "short_candidate": parent["candidate"],
                   "short_loss": parent["loss"], **record, **check}
            deep_rows.append(row)
            checkpoint(short_rows,deep_rows)
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
              "recovery": "seven run351 gate lists reused, optimizer metadata unavailable; 53 new short fits; persistent checkpoint permits resume without overwriting originals",
              "short_count": len(short_rows), "deep_count": len(deep_rows),
              "short_rows": short_rows, "deep_rows": deep_rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "sixty specified one-slot mutations of opening02 source, each making two changes from original chain; one warm start, six selected by fitted loss; bounded numerical evidence only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
