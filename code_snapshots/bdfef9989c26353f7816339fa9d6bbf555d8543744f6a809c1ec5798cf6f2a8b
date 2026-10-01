"""Closing-edge mutations of the fresh run339 chain13 numerical lead.

Replace each of the thirteen CNOT slots by 0->3 at even slots and 3->0
at odd slots. The source uses only 01,12,23 chain edges. All 168 local
angles remain free; there are no inherited fixed local chunks.
"""

import argparse
import hashlib
import json
from pathlib import Path

import torch

from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_prefix_perturb import decode_layers


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "search13_fresh_chain_refine_trial1.json"
PRIOR = ROOT / "search13_fresh_chain_refine_result.json"
OUT = ROOT / "search13_fresh_chain_closing_edge_result.json"
torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=300)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    angles = decode_layers(source)
    schedule = tuple((g["control"], g["target"])
                     for g in source["gates"] if g["gate"] == "cx")
    assert len(schedule) == 13
    assert {tuple(sorted(edge)) for edge in schedule} == {(0, 1), (1, 2), (2, 3)}
    assert len(source["gates"]) == 14 * 4 + 13
    assert all(g["gate"] in ("u3", "cx") for g in source["gates"])
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    chunks = [[] for _ in range(14)]
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed,
                          matrix_for_schedule(schedule), v, False)
    expected = next(row["loss"] for row in json.loads(PRIOR.read_text())["rows"]
                    if row["trial"] == 1)
    assert abs(roundtrip - expected) < 1e-10, (roundtrip, expected)
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "expected_loss": expected,
                      "roundtrip_error": abs(roundtrip - expected)}), flush=True)
    rows = []
    for slot in range(13):
        mutated = list(schedule)
        replacement = (0, 3) if slot % 2 == 0 else (3, 0)
        mutated[slot] = replacement
        fitted, record = fit(mutated, angles, fixed, v, args.maxiter)
        candidate = serialize(chunks, mutated, fitted)
        path = ROOT / f"search13_fresh_chain_closing_edge_slot{slot}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"slot": slot, "old": schedule[slot],
               "replacement": replacement, "schedule": mutated, **record,
               "candidate": path.name, "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps({"slot": slot, "old": row["old"],
                          "replacement": replacement,
                          "loss": row["loss"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
    assert len({tuple(map(tuple, row["schedule"])) for row in rows}) == 13
    best = min(rows, key=lambda row: row["off_diagonal_error"])
    result = {"source": SOURCE.name,
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "source_schedule": schedule, "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "fixed_local_matrices": "identity at all 14 layers; no inherited chunks",
              "proposal_rule": "replace each slot by closing edge 03, 0->3 at even slots and 3->0 at odd slots",
              "maxiter": args.maxiter, "fits": len(rows), "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "thirteen specified one-slot closing-edge changes, one source warm start each, all 168 local parameters free; failure is not a no-go"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
