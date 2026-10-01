"""Joint fit for endpoint-guided five-CNOT prefixes and eight-CNOT tail.

The selected prefix endpoints differ by one CNOT on a first matchgate pair
from the exact six-CNOT endpoint.  They cover two exact14 nonlocal phase
masks.  All local SU(2) gates before, at, and after the boundary are free;
the exact eight-CNOT tail edge pattern supplies only an interaction schedule.
"""

import argparse
import json

import numpy as np

from check_circuit import evaluate, one_qubit_matrix
from search13_joint_gauge_five_endpoint import (
    ROOT, WARM, candidate, endpoint, fit, local_to_axis, warm_angles,
)


SCREEN = ROOT / "search13_endpoint_guided_screen_result.json"
BASE = ROOT / "topology14_exact_matchgate_rational.json"
SELECTED = (("2->0", 0), ("2->0", 3), ("1->3", 0), ("1->3", 1))
MASK_PHASE = {7: -np.pi/4, 6: -np.pi/2,
              14: np.pi/8, 11: -np.pi/8}


def tail_chunks_and_edges():
    gates = json.loads(BASE.read_text())["gates"][13:]
    chunks, edges = [[]], []
    for gate in gates:
        if gate["gate"] == "cx":
            edges.append((gate["control"], gate["target"]))
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 9 and len(edges) == 8
    return chunks, tuple(edges)


def phase_informed_angles(prefix, tail_chunks):
    mats = [[np.eye(2, dtype=complex) for _ in range(4)] for _ in range(14)]

    def apply(slot, gate):
        q = gate["qubit"]
        mats[slot][q] = one_qubit_matrix(gate) @ mats[slot][q]

    for q, theta in ((1, "5*pi/8"), (2, "3*pi/4"), (3, "3*pi/8")):
        apply(0, {"gate": "rz", "qubit": q, "theta": theta})
    rows = [8, 4, 2, 1]
    applied = set()
    for index, (a, b) in enumerate(prefix):
        rows[b] ^= rows[a]
        mask = rows[b]
        if mask in MASK_PHASE and mask not in applied:
            apply(index+1, {"gate": "rz", "qubit": b,
                            "theta": MASK_PHASE[mask]})
            applied.add(mask)
    for index, chunk in enumerate(tail_chunks):
        for gate in chunk:
            apply(index+5, gate)
    return np.array([[local_to_axis(mats[slot][wire]) for wire in range(4)]
                     for slot in range(14)]), tuple(sorted(applied))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=25050)
    parser.add_argument("--maxiter", type=int, default=350)
    parser.add_argument("--phase-sigma", type=float, default=0.04)
    parser.add_argument("--warm-sigma", type=float, default=0.08)
    args = parser.parse_args()
    screen = json.loads(SCREEN.read_text())
    tail_chunks, tail_edges = tail_chunks_and_edges()
    warm = warm_angles()
    rng = np.random.default_rng(args.seed)
    rows = []
    for index, (orientation, rank) in enumerate(SELECTED):
        entry = screen["selected"][orientation][rank]
        prefix = tuple(tuple(edge) for edge in entry["path"])
        assert endpoint(prefix) == tuple(entry["endpoint"])
        seeded, applied = phase_informed_angles(prefix, tail_chunks)
        assert len(applied) == 2
        edges = prefix+tail_edges
        for kind, start, sigma in (("phase", seeded, args.phase_sigma),
                                   ("warm", warm, args.warm_sigma)):
            initial = start+rng.normal(0, sigma, start.shape)
            angles, row = fit(edges, initial, args.maxiter)
            circuit = candidate(edges, angles)
            path = ROOT / f"search13_endpoint_guided_case{index}_{kind}.json"
            path.write_text(json.dumps(circuit, indent=2)+"\n")
            checked = evaluate(circuit)
            row.update(case=index, start=kind, sigma=sigma,
                       orientation=orientation, prefix=prefix,
                       endpoint=entry["endpoint"],
                       visited_phase_masks=applied,
                       omitted_phase_masks=entry["omitted_required_masks"],
                       tail_edges=tail_edges, candidate=path.name,
                       off_diagonal_error=checked["off_diagonal_error"],
                       valid_diagonalizer=checked["valid_diagonalizer"])
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    result = {"screen": SCREEN.name, "source": BASE.name,
              "warm_source": WARM.name, "selected": SELECTED,
              "seed": args.seed, "maxiter": args.maxiter,
              "phase_sigma": args.phase_sigma, "warm_sigma": args.warm_sigma,
              "rows": rows,
              "best_loss": min(rows, key=lambda r: r["loss"]),
              "best_off_diagonal": min(rows, key=lambda r: r["off_diagonal_error"])}
    (ROOT / "search13_endpoint_guided_fit_result.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
