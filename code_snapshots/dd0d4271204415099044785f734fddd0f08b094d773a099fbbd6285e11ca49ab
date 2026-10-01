"""Six diverse seven-CNOT tail fits for the refined numerical output gauge.

Select one schedule from each distance-to-one-deletion class 2..7 among the
numerically rank-compatible schedules. Fit all local SU(2) layers first to
the fixed gauged tail by process overlap, then directly to the label-free
cycle diagonalization condition with the exact six-CNOT prefix fixed.
"""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import torch

from check_circuit import evaluate
from search13_eigenrow_gauge_seven_fit import fit as process_fit
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import gauge, setup
from search13_fulltail_rank8_topology import edit_distance_to_one_deletion


ROOT = Path(__file__).resolve().parent
TOPOLOGY = ROOT / "search13_refined_gauge_topology_result.json"
GAUGE = ROOT / "search13_mixedcut_gauge_refine_result.json"
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_refined_gauge_fit_result.json"
torch.set_num_threads(1)


def choose_schedules(survivors):
    groups = {d: [] for d in range(2, 8)}
    for raw in survivors:
        schedule = tuple(tuple(pair) for pair in raw)
        distance = edit_distance_to_one_deletion(schedule)
        if distance in groups:
            groups[distance].append(schedule)
    selected = []
    for distance in range(2, 8):
        assert groups[distance]

        def score(schedule):
            counts = Counter(schedule)
            nearest = min((sum(a != b for a, b in zip(schedule, old))
                           for old in selected), default=7)
            adjacent_repeats = sum(schedule[j] == schedule[j+1]
                                   for j in range(6))
            return (nearest, len(counts), -max(counts.values()),
                    -adjacent_repeats, schedule)

        selected.append(max(groups[distance], key=score))
    return selected, {str(d): len(rows) for d, rows in groups.items()}


def directed(schedule, case):
    # Direction reversal is locally equivalent but offers distinct starts.
    return tuple(pair[::-1] if (case + slot) % 2 else pair
                 for slot, pair in enumerate(schedule))


def write_checked(name, prefix, tail_gates):
    candidate = {"n": 4, "gates": prefix + tail_gates}
    path = ROOT / name
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    checked = evaluate(candidate)
    return {"candidate": path.name,
            "cnot_count": checked["cnot_count"],
            "off_diagonal_error": checked["off_diagonal_error"],
            "valid_diagonalizer": checked["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=30100)
    parser.add_argument("--process-maxiter", type=int, default=200)
    parser.add_argument("--direct-maxiter", type=int, default=300)
    args = parser.parse_args()
    topology = json.loads(TOPOLOGY.read_text())
    gauge_data = json.loads(GAUGE.read_text())
    selected, counts = choose_schedules(topology["surviving_schedules"])
    tail, prefix_matrix, labels, sectors = setup()
    assert len(gauge_data["refined_angles"]) == sum(len(s)**2 for s in sectors)
    h = gauge(torch.as_tensor(gauge_data["refined_angles"], dtype=torch.float64),
              sectors).detach().numpy()
    assert np.max(np.abs(h @ h.conj().T - np.eye(16))) < 1e-10
    assert all(abs(h[i, j]) < 1e-10 or labels[i] == labels[j]
               for i in range(16) for j in range(16))
    target = h @ tail
    prefix_gates = json.loads(SOURCE.read_text())["gates"][:13]
    rows = []
    for case, unoriented in enumerate(selected):
        edges = directed(unoriented, case)
        seed = args.seed + 137*case
        process_gates, process = process_fit(edges, seed,
                                              args.process_maxiter, target)
        process.update(write_checked(
            f"search13_refined_gauge_fit_case{case}_process.json",
            prefix_gates, process_gates))
        angles = decode_local_layers(process_gates)
        direct_gates, direct = direct_fit(edges, prefix_matrix, angles,
                                         args.direct_maxiter)
        direct.update(write_checked(
            f"search13_refined_gauge_fit_case{case}_direct.json",
            prefix_gates, direct_gates))
        row = {"case": case,
               "distance_to_one_deletion": edit_distance_to_one_deletion(unoriented),
               "unoriented_schedule": unoriented, "directed_edges": edges,
               "seed": seed, "process": process, "direct": direct}
        rows.append(row)
        print(json.dumps({"case": case,
                          "distance": row["distance_to_one_deletion"],
                          "process_loss": process["loss"],
                          "direct_loss": direct["loss"],
                          "direct_offdiag": direct["off_diagonal_error"],
                          "direct_valid": direct["valid_diagonalizer"]},
                         sort_keys=True), flush=True)
    result = {"source": SOURCE.name, "topology_source": TOPOLOGY.name,
              "gauge_source": GAUGE.name,
              "numerical_rank_threshold": topology["numerical_rank_threshold"],
              "survivor_distance_counts": counts, "selected_schedules": selected,
              "seed": args.seed, "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter, "rows": rows,
              "best_direct_case": min(rows, key=lambda r: r["direct"]["loss"])["case"],
              "scope": "six rank-compatible topologies for an approximate output gauge; numerical fits only"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
