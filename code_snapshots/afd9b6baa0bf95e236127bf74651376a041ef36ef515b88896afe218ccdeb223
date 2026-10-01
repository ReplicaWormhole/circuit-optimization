"""Alternating exact block-Procrustes gauge and seven-CNOT local-layer fits.

The first six CNOTs of the exact 14-CNOT circuit are fixed. For a selected
seven-CNOT tail schedule, each outer iteration analytically chooses the best
output unitary in the entire commuting eigenspace gauge U(6)xU(4)xU(3)xU(3),
then optimizes the tail's 96 axis-angle local parameters for that fixed
gauged target. This is a bounded nonconvex numerical search.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from check_circuit import evaluate
from search13_fulltail_invariant_fit import decode_local_layers, direct_fit
from search13_fulltail_invariant_gauge_rank import circuit_matrix, setup
from search13_identity_tail_fit import (
    choose, deletion_warm, original_tail, process_fit,
)


ROOT = Path(__file__).resolve().parent
TOPOLOGY = ROOT / "search13_identity_tail_topology_result.json"
OUT = ROOT / "search13_procrustes_gauge_fit_result.json"


def optimal_gauge(unitary, target, sectors):
    """Return the block unitary minimizing ||unitary-G target||_F."""
    cross = unitary @ target.conj().T
    gauge = np.zeros((16, 16), dtype=complex)
    for sector in sectors:
        indices = np.ix_(sector, sector)
        left, _, right = np.linalg.svd(cross[indices], full_matrices=False)
        gauge[indices] = left @ right
    return gauge


def orbit_loss(unitary, target, gauge):
    return float(np.linalg.norm(unitary - gauge @ target) ** 2 / 16)


def save_candidate(path, prefix, tail_gates):
    candidate = {"n": 4, "gates": prefix + tail_gates}
    path.write_text(json.dumps(candidate, indent=2) + "\n")
    check = evaluate(candidate)
    return {"candidate": path.name, "cnot_count": check["cnot_count"],
            "off_diagonal_error": check["off_diagonal_error"],
            "valid_diagonalizer": check["valid_diagonalizer"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=32600)
    parser.add_argument("--outer", type=int, default=3)
    parser.add_argument("--process-maxiter", type=int, default=60)
    parser.add_argument("--direct-maxiter", type=int, default=120)
    args = parser.parse_args()

    topology = json.loads(TOPOLOGY.read_text())
    selected, _ = choose(topology["surviving_schedules"], 6)
    cases = (0, 5)
    prefix, original_gates, original_pairs = original_tail()
    tail, prefix_matrix, labels, sectors = setup()
    rows = []
    for case in cases:
        edges = selected[case]
        warm, deleted, distance = deletion_warm(
            original_gates, original_pairs, edges)
        rng = np.random.default_rng(args.seed + case)
        starts = (("deletion_warm", warm),
                  ("random", rng.normal(0, 0.4, (8, 4, 3))))
        for start_name, initial in starts:
            # The first matrix is a genuine seven-CNOT circuit with local
            # rotations at the deletion warm start; process_fit serializes it.
            current_gates, _ = process_fit(edges, initial, tail, 0)
            stages = []
            for outer in range(args.outer):
                before = circuit_matrix(current_gates)
                gauge = optimal_gauge(before, tail, sectors)
                assert np.max(np.abs(gauge @ gauge.conj().T - np.eye(16))) < 1e-12
                assert all(abs(gauge[i, j]) < 1e-12 or labels[i] == labels[j]
                           for i in range(16) for j in range(16))
                before_loss = orbit_loss(before, tail, gauge)
                current_gates, process = process_fit(
                    edges, decode_local_layers(current_gates),
                    gauge @ tail, args.process_maxiter)
                after = circuit_matrix(current_gates)
                next_gauge = optimal_gauge(after, tail, sectors)
                stages.append({"outer": outer, "process": process,
                               "orbit_loss_before": before_loss,
                               "orbit_loss_after": orbit_loss(after, tail, next_gauge)})
            process_path = ROOT / f"search13_procrustes_case{case}_{start_name}_process.json"
            process_check = save_candidate(process_path, prefix, current_gates)
            direct_gates, direct = direct_fit(
                edges, prefix_matrix, decode_local_layers(current_gates),
                args.direct_maxiter)
            direct_path = ROOT / f"search13_procrustes_case{case}_{start_name}_direct.json"
            direct_check = save_candidate(direct_path, prefix, direct_gates)
            row = {"case": case, "schedule": edges, "start": start_name,
                   "deleted_exact_tail_cnot": deleted,
                   "distance_to_deletion": distance, "stages": stages,
                   "process_check": process_check,
                   "direct": {**direct, **direct_check}}
            rows.append(row)
            print(json.dumps({"case": case, "start": start_name,
                              "orbit_loss": stages[-1]["orbit_loss_after"],
                              "direct_loss": direct["loss"],
                              "off_diagonal_error": direct_check["off_diagonal_error"],
                              "valid": direct_check["valid_diagonalizer"]}),
                  flush=True)
    best = min(rows, key=lambda r: r["direct"]["off_diagonal_error"])
    result = {"source": "topology14_exact_matchgate_rational.json",
              "topology_source": TOPOLOGY.name,
              "strategy": "alternating block-Procrustes output eigenspace gauge and fixed-target process fit, followed by label-free direct fit",
              "gauge_sector_dimensions": list(map(len, sectors)),
              "seed": args.seed, "outer": args.outer,
              "process_maxiter": args.process_maxiter,
              "direct_maxiter": args.direct_maxiter,
              "case_indices": cases, "rows": rows,
              "best_candidate": best["direct"]["candidate"],
              "best_off_diagonal_error": best["direct"]["off_diagonal_error"],
              "scope": "four bounded numerical fits with the exact six-CNOT prefix fixed; failure does not exclude other gauges or circuits"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
