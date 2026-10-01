"""Bounded numerical fits of six diverse distance-three rank-eight tail orders.

The selected orders survived the fixed-gauge all-cut exact necessary screen.
For each, fit the exact gauged tail by process overlap, continue using the
label-free V4 objective, and independently fit that objective from a structured
one-deletion/three-rewire exact-eight-tail warm start. Numerical evidence only.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate, one_qubit_matrix
from search13_eigenrow_gauge_seven_fit import fit as process_fit
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import setup
from search13_fulltail_rank8_fit import exact_gauge_numeric
from search13_fulltail_rank8_neighborhood_fit import split_tail
from search13_joint_gauge_five_endpoint import local_to_axis


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
SCREEN = ROOT / "search13_rank8_global_allcuts_result.json"
OUT = ROOT / "search13_fulltail_rank8_distance3_fit_result.json"
SELECTED_INDICES = (341, 0, 295, 446, 248, 1459)
CUTS = ({0, 1}, {0, 2}, {0, 3})


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def profile(schedule):
    return tuple(sum((a in cut) != (b in cut) for a, b in schedule)
                 for cut in CUTS)


def warm_start(schedule):
    chunks, old_cxs = split_tail()
    options = []
    for deleted in range(8):
        retained = [edge for i, edge in enumerate(old_cxs) if i != deleted]
        distance = hamming([tuple(sorted(edge)) for edge in retained],
                           [tuple(edge) for edge in schedule])
        options.append((distance, deleted, retained))
    distance, deleted, retained = min(options, key=lambda x: x[:2])
    assert distance == 3
    mismatch = [i for i, (old, new) in enumerate(zip(retained, schedule))
                if tuple(sorted(old)) != tuple(new)]
    assert len(mismatch) == 3
    directed = tuple(old if i not in mismatch else tuple(new)
                     for i, (old, new) in enumerate(zip(retained, schedule)))
    merged = [list(chunk) for chunk in chunks]
    merged[deleted].extend(merged[deleted + 1])
    del merged[deleted + 1]
    assert len(merged) == 8
    angles = np.zeros((8, 4, 3))
    for slot, chunk in enumerate(merged):
        for qubit in range(4):
            local = np.eye(2, dtype=complex)
            for gate in chunk:
                if gate["qubit"] == qubit:
                    local = one_qubit_matrix(gate) @ local
            angles[slot, qubit] = local_to_axis(local)
    return directed, angles, {"deleted_old_tail_cx": deleted,
                              "rewired_retained_indices": mismatch}


def save(prefix, tail, case, stage):
    candidate = {"n": 4, "gates": prefix + tail}
    path = ROOT / f"search13_fulltail_rank8_distance3_fit_case{case}_{stage}.json"
    assert not path.exists(), f"refusing overwrite of {path}"
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process-maxiter", type=int, default=220)
    parser.add_argument("--direct-maxiter", type=int, default=320)
    parser.add_argument("--seed", type=int, default=29371)
    args = parser.parse_args()
    assert not OUT.exists(), f"refusing overwrite of {OUT}"
    source = json.loads(SCREEN.read_text())
    schedules = [tuple(map(tuple, source["surviving_schedules"][i]))
                 for i in SELECTED_INDICES]
    assert len(schedules) == len(set(schedules)) == 6
    assert {profile(s) for s in schedules} == {(4, 4, 6), (4, 5, 5),
                                               (5, 4, 5), (5, 5, 4)}
    assert all(all(s[i] != s[i + 1] for i in range(6)) for s in schedules)
    assert min(hamming(x, y) for i, x in enumerate(schedules)
               for y in schedules[i + 1:]) >= 5
    tail, prefix_matrix, _, _ = setup()
    target = exact_gauge_numeric() @ tail
    prefix = json.loads(SOURCE.read_text())["gates"][:13]
    rows = []
    for case, (index, schedule) in enumerate(zip(SELECTED_INDICES, schedules)):
        directed, warm_angles, alignment = warm_start(schedule)
        process_gates, process = process_fit(directed, args.seed + case,
                                             args.process_maxiter, target)
        process.update(save(prefix, process_gates, case, "process"))
        direct_gates, direct = direct_fit(directed, prefix_matrix,
                                          decode_local_layers(process_gates),
                                          args.direct_maxiter)
        direct.update(save(prefix, direct_gates, case, "direct"))
        warm_gates, warm = direct_fit(directed, prefix_matrix, warm_angles,
                                      args.direct_maxiter)
        warm.update(save(prefix, warm_gates, case, "warm_direct"))
        row = {"case": case, "survivor_index": index, "schedule": schedule,
               "directed_schedule": directed, "balanced_crossings": profile(schedule),
               "alignment": alignment, "process": process, "direct": direct,
               "warm_direct": warm}
        rows.append(row)
        print(json.dumps({"case": case, "profile": profile(schedule),
                          "process_loss": process["loss"],
                          "direct_loss": direct["loss"], "warm_loss": warm["loss"],
                          "best_offdiag": min(direct["off_diagonal_error"],
                                              warm["off_diagonal_error"])}),
              flush=True)
    OUT.write_text(json.dumps({"source": SOURCE.name, "screen": SCREEN.name,
                               "exact_gauge": "search13_fulltail_invariant_exact_witness_result.json",
                               "selected_indices": SELECTED_INDICES,
                               "selection": "cover all four balanced crossing profiles, zero adjacent pair repeats, pairwise Hamming distance >=5",
                               "seed": args.seed, "process_maxiter": args.process_maxiter,
                               "direct_maxiter": args.direct_maxiter,
                               "rows": rows,
                               "scope": "six fixed seven-CNOT tail orders and one direction assignment each; one random process start and two direct starts; numerical screen, not a no-go"},
                              indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
