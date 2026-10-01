"""Bounded numerical fits for the eight mixed-cut-compatible tail schedules.

Each schedule receives a random-start process fit to the exact rank-eight
gauged tail, a label-free continuation, and a second label-free fit seeded by
deleting/rewiring the exact eight-CNOT tail. This is a numerical screen only.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate, one_qubit_matrix
from delete14_search import fixed_matrix
from search13_eigenrow_gauge_seven_fit import fit as process_fit
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_rank8_fit import exact_gauge_numeric
from search13_fulltail_invariant_gauge_rank import setup
from search13_joint_gauge_five_endpoint import local_to_axis


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
SCREEN = ROOT / "search13_rank8_neighborhood_mixedcut_result.json"
OUT = ROOT / "search13_fulltail_rank8_neighborhood_fit_result.json"


def split_tail():
    chunks = [[]]
    cxs = []
    for gate in json.loads(SOURCE.read_text())["gates"][13:]:
        if gate["gate"] == "cx":
            cxs.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(cxs) == 8 and len(chunks) == 9
    return chunks, cxs


def warm_start(schedule):
    """Choose an exact-eight-tail deletion with at most one pair mismatch."""
    chunks, cxs = split_tail()
    normalized = lambda edge: tuple(sorted(edge))
    options = []
    for deleted in range(8):
        retained = [edge for i, edge in enumerate(cxs) if i != deleted]
        mismatch = [i for i, (old, new) in enumerate(zip(retained, schedule))
                    if normalized(old) != normalized(new)]
        if len(mismatch) <= 1:
            options.append((len(mismatch), deleted, retained, mismatch))
    assert options, f"schedule is outside one-deletion/one-rewire neighborhood: {schedule}"
    _, deleted, retained, mismatch = min(options, key=lambda item: item[:2])
    directed = [old if i not in mismatch else tuple(new)
                for i, (old, new) in enumerate(zip(retained, schedule))]
    merged = [list(chunk) for chunk in chunks]
    merged[deleted].extend(merged[deleted + 1])
    del merged[deleted + 1]
    assert len(merged) == 8 and len(directed) == 7
    angles = np.zeros((8, 4, 3))
    for slot, chunk in enumerate(merged):
        for qubit in range(4):
            local = np.eye(2, dtype=complex)
            for gate in chunk:
                if gate["qubit"] == qubit:
                    local = one_qubit_matrix(gate) @ local
            angles[slot, qubit] = local_to_axis(local)
    return tuple(directed), angles, {"deleted_old_tail_cx": deleted,
                                     "rewired_retained_index": mismatch[0] if mismatch else None}


def store(prefix_gates, tail_gates, case, stage):
    candidate = {"n": 4, "gates": prefix_gates + tail_gates}
    path = ROOT / f"search13_fulltail_rank8_neighborhood_fit_case{case}_{stage}.json"
    assert not path.exists(), f"refusing to overwrite {path}"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name,
            "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process-maxiter", type=int, default=140)
    parser.add_argument("--direct-maxiter", type=int, default=220)
    parser.add_argument("--seed", type=int, default=28500)
    args = parser.parse_args()
    assert not OUT.exists(), f"refusing to overwrite {OUT}"
    schedules = json.loads(SCREEN.read_text())["surviving_schedules"]
    assert len(schedules) == 8
    tail, prefix_matrix, _, _ = setup()
    target = exact_gauge_numeric() @ tail
    prefix_gates = json.loads(SOURCE.read_text())["gates"][:13]
    rows = []
    for case, schedule in enumerate(schedules):
        directed, warm, alignment = warm_start(schedule)
        process_gates, process = process_fit(directed, args.seed + case,
                                             args.process_maxiter, target)
        process.update(store(prefix_gates, process_gates, case, "process"))
        process_angles = decode_local_layers(process_gates)
        direct_gates, direct = direct_fit(directed, prefix_matrix,
                                          process_angles, args.direct_maxiter)
        direct.update(store(prefix_gates, direct_gates, case, "direct"))
        warm_gates, warm_direct = direct_fit(directed, prefix_matrix,
                                             warm, args.direct_maxiter)
        warm_direct.update(store(prefix_gates, warm_gates, case, "warm_direct"))
        row = {"case": case, "undirected_schedule": schedule,
               "directed_schedule": directed, "alignment": alignment,
               "process": process, "direct": direct,
               "warm_direct": warm_direct}
        rows.append(row)
        print(json.dumps({"case": case, "process_loss": process["loss"],
                          "direct_loss": direct["loss"],
                          "warm_direct_loss": warm_direct["loss"],
                          "best_offdiag": min(direct["off_diagonal_error"],
                                              warm_direct["off_diagonal_error"]) }),
              flush=True)
    OUT.write_text(json.dumps({"source": SOURCE.name, "screen": SCREEN.name,
                               "exact_gauge": "search13_fulltail_invariant_exact_witness_result.json",
                               "seed": args.seed, "process_maxiter": args.process_maxiter,
                               "direct_maxiter": args.direct_maxiter, "rows": rows,
                               "scope": "eight specified physical-pair orders; inherited directions for unchanged gates, ascending directions for rewired gates; one random process plus two direct starts each; numerical only"},
                              indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
