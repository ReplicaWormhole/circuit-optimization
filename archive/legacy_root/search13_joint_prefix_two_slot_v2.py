"""Twenty deterministic coordinated two-slot mutations from run332's basin.

Proposals cover ten prefix/tail slot pairs, seven adjacent tail pairs, and
three separated tail pairs. Prefer disjoint wire pairs, then use a
deterministic support fallback to exclude prior unordered schedules.
Directions are mixed. No proposal is selected by raw loss.
All 168 local angles are fitted jointly for at most 300 iterations.
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
SOURCE = ROOT / "search13_joint_prefix_slot8_refine_trial1.json"
PRIOR = ROOT / "search13_joint_prefix_slot8_refine_result.json"
OUT = ROOT / "search13_joint_prefix_two_slot_v2_result.json"
SLOT_PAIRS = (
    (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
    (0, 10), (1, 11), (2, 12), (3, 6), (4, 7),
    (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12),
    (5, 9), (6, 10), (7, 11),
)
torch.set_num_threads(1)


def topology(candidate):
    return tuple((g["control"], g["target"]) for g in candidate["gates"]
                 if g["gate"] == "cx")


def prior_topologies():
    out = set()
    for path in ROOT.glob("search13*.json"):
        data = json.loads(path.read_text())
        if isinstance(data, dict) and data.get("n") == 4 and "gates" in data:
            out.add(topology(data))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=300)
    args = parser.parse_args()
    assert not OUT.exists(), OUT
    source = json.loads(SOURCE.read_text())
    angles = decode_layers(source)
    expected = next(row["loss"] for row in json.loads(PRIOR.read_text())["rows"]
                    if row["trial"] == 1)
    deleted, rewire, blocks = CASES[1]
    chunks, cxs, fixed, _ = prepare(deleted, blocks, prefix_rewire=rewire)
    schedule = [(g["control"], g["target"]) for g in cxs]
    schedule[8] = (1, 2)
    assert tuple(schedule) == topology(source)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    roundtrip = objective(angles.ravel(), fixed,
                          matrix_for_schedule(schedule), v, False)
    assert abs(roundtrip - expected) < 1e-10, (roundtrip, expected)
    print(json.dumps({"source_roundtrip_loss": roundtrip,
                      "source_expected_loss": expected,
                      "roundtrip_error": abs(roundtrip - expected)}), flush=True)
    excluded = prior_topologies()
    excluded_unordered = {tuple(tuple(sorted(edge)) for edge in s)
                          for s in excluded}
    proposals = []
    for case, slots in enumerate(SLOT_PAIRS):
        options = []
        for offset, slot in enumerate(slots):
            old = schedule[slot]
            complement = tuple(q for q in range(4) if q not in old)
            others = [edge for edge in ((a, b) for a in range(4)
                                        for b in range(a+1, 4))
                      if set(edge) != set(old) and edge != complement]
            # Prefer disjoint replacement. Fall back to a different support
            # only if that pair of replacements reproduces an old schedule.
            base = [complement] + others
            options.append([edge[::-1] if (case + offset) % 2 else edge
                            for edge in base])
        chosen = None
        for left in options[0]:
            for right in options[1]:
                mutated = list(schedule)
                mutated[slots[0]], mutated[slots[1]] = left, right
                unordered = tuple(tuple(sorted(edge)) for edge in mutated)
                if (tuple(mutated) not in excluded and
                        unordered not in excluded_unordered and
                        all(tuple(mutated) != tuple(s) for s, _ in proposals)):
                    chosen = mutated, [{"slot": slot, "old": schedule[slot],
                                        "replacement": edge}
                                       for slot, edge in zip(slots, (left, right))]
                    break
            if chosen is not None:
                break
        assert chosen is not None, (case, slots)
        mutated, changes = chosen
        assert sum(a != b for a, b in zip(mutated, schedule)) == 2
        assert tuple(mutated) not in excluded, (case, mutated)
        proposals.append((mutated, changes))
    assert len({tuple(s) for s, _ in proposals}) == 20
    rows = []
    for case, (mutated, changes) in enumerate(proposals):
        fitted, record = fit(mutated, angles, fixed, v, args.maxiter)
        candidate = serialize(chunks, mutated, fitted)
        path = ROOT / f"search13_joint_prefix_two_slot_v2_case{case}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        row = {"case": case, "changes": changes, **record,
               "candidate": path.name, "cnot_count": check["cnot_count"],
               "off_diagonal_error": check["off_diagonal_error"],
               "valid_diagonalizer": check["valid_diagonalizer"]}
        rows.append(row)
        print(json.dumps({"case": case, "slots": SLOT_PAIRS[case],
                          "loss": row["loss"],
                          "off_diagonal_error": row["off_diagonal_error"],
                          "valid": row["valid_diagonalizer"]}), flush=True)
    best = min(rows, key=lambda row: row["off_diagonal_error"])
    result = {"source": SOURCE.name, "source_roundtrip_loss": roundtrip,
              "source_roundtrip_error": abs(roundtrip - expected),
              "source_topology": schedule, "slot_pairs": SLOT_PAIRS,
              "excluded_prior_directed_topologies": len(excluded),
              "excluded_prior_unordered_topologies": len(excluded_unordered),
              "proposal_rule": "prefer complementary disjoint replacements, with deterministic support fallback to exclude both previously fitted directed and unordered schedules; mixed orientations; no raw-loss selection",
              "fits": len(rows), "maxiter": args.maxiter, "rows": rows,
              "best_candidate": best["candidate"],
              "best_off_diagonal_error": best["off_diagonal_error"],
              "scope": "twenty specified two-slot topologies, one fixed warm start each, full local freedom; numerical failures do not exclude topologies"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
