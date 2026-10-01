"""Fit diverse one-slot topology changes from the refined run-327 near-hit.

For each of all thirteen CNOT slots, replace its interaction by the disjoint
wire pair; alternate orientations by slot parity. Also reverse directions
at slots 0,4,8,12. Proposal selection is independent of numerical raw loss.
Every fit varies all 168 local SU(2) parameters from the same source warm
start. These are bounded numerical fits, not lower-bound evidence.
"""

import argparse
import json
from pathlib import Path

import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_joint_prefix_orbit_fit import CASES
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_joint_prefix_orbit_refine_0.json"
PRIOR = ROOT / "search13_joint_prefix_one_step_escape_result.json"
OUT = ROOT / "search13_joint_prefix_slot_diversity_result.json"
torch.set_num_threads(1)


def proposals(schedule):
    for slot, old in enumerate(schedule):
        replacement = tuple(q for q in range(4) if q not in old)
        if slot % 2:
            replacement = replacement[::-1]
        yield {"kind": "disjoint_pair", "slot": slot, "old": old,
               "replacement": replacement}
    for slot in (0, 4, 8, 12):
        old = schedule[slot]
        yield {"kind": "direction_reversal", "slot": slot, "old": old,
               "replacement": old[::-1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=250)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    angles = decode_layers(source)
    deleted, rewire, blocks = CASES[1]
    chunks, old_cxs, fixed, _ = prepare(deleted, blocks,
                                           prefix_rewire=rewire)
    schedule = tuple((g["control"], g["target"]) for g in old_cxs)
    assert tuple((g["control"], g["target"])
                 for g in source["gates"] if g["gate"] == "cx") == schedule
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    baseline = objective(angles.ravel(), fixed,
                         matrix_for_schedule(schedule), v, False)
    assert abs(baseline - 0.0787330125161619) < 1e-8, baseline
    selected = list(proposals(schedule))
    assert len(selected) == 17
    assert {r["slot"] for r in selected} == set(range(13))
    prior = json.loads(PRIOR.read_text())
    previously_fitted = {(r["slot"], tuple(r["replacement"]))
                         for r in prior["rows"]}
    assert not any((r["slot"], r["replacement"]) in previously_fitted
                   for r in selected)
    rows = []
    for case, proposal in enumerate(selected):
        mutated = list(schedule)
        mutated[proposal["slot"]] = proposal["replacement"]
        raw = objective(angles.ravel(), fixed,
                        matrix_for_schedule(mutated), v, False)
        fitted, record = fit(mutated, angles, fixed, v, args.maxiter)
        candidate = serialize(chunks, mutated, fitted)
        path = ROOT / f"search13_joint_prefix_slot_diversity_case{case}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"case": case, **proposal, "raw_loss": raw, **record,
               "candidate": path.name, "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps({"case": case, "kind": row["kind"],
                          "slot": row["slot"], "old": row["old"],
                          "replacement": row["replacement"],
                          "loss": row["loss"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["off_diagonal_error"])
    result = {"source": SOURCE.name, "source_topology": schedule,
              "baseline_loss": baseline,
              "selection_rule": "disjoint wire pair at every slot with alternating orientation, plus direction reversals at slots 0,4,8,12; no raw-loss ranking",
              "maxiter": args.maxiter, "fits": len(rows), "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "17 one-step directed topology fits from one fixed near-hit warm start; all 168 local layers free; failures are not exclusions"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
