"""Bounded seven-CNOT fixed-unitary fit for a pi/4 eigenrow Givens gauge.

Rows 8 and 11 have the same V4 eigenvalue. Their real pi/4 Givens rotation
preserves exact diagonalization and gives full-tail balanced ranks
(16,13,16). Fit the two nearest rank-compatible seven-CNOT schedules used
for the row-swap comparison, one fixed seed each.
"""

import argparse
import json

import numpy as np

from check_circuit import evaluate
from delete14_search import fixed_matrix
from search13_eigenrow_gauge_seven_fit import BASE, ROOT, SCHEDULES, fit


def target_and_prefix():
    gates = json.loads(BASE.read_text())["gates"]
    prefix, tail = gates[:13], gates[13:]
    target = np.eye(16, dtype=complex)
    for gate in tail:
        target = fixed_matrix(gate) @ target
    a, b = target[8].copy(), target[11].copy()
    co = np.cos(np.pi / 4)
    si = np.sin(np.pi / 4)
    target[8] = co*a - si*b
    target[11] = si*a + co*b
    return prefix, target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxiter", type=int, default=250)
    args = parser.parse_args()
    prefix, target = target_and_prefix()
    rows = []
    for schedule_index, edges in enumerate(SCHEDULES):
        seed = 2301 + schedule_index
        tail, record = fit(edges, seed, args.maxiter, target)
        candidate = {"n": 4, "gates": list(prefix) + tail}
        path = ROOT / f"search13_eigenrow_givens_seven_fit_top{schedule_index}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        check = evaluate(candidate)
        record.update(schedule_index=schedule_index, edges=edges,
                      candidate=path.name,
                      off_diagonal_error=check["off_diagonal_error"],
                      valid_diagonalizer=check["valid_diagonalizer"])
        rows.append(record)
        print(json.dumps(record, sort_keys=True), flush=True)
    result = {"gauge": {"kind": "real-givens", "rows": [8, 11],
                        "angle_pi_multiple": "1/4"},
              "schedules": SCHEDULES, "maxiter": args.maxiter,
              "rows": rows,
              "best_process_loss": min(rows, key=lambda row: row["loss"]),
              "best_cycle_error": min(rows, key=lambda row: row["off_diagonal_error"])}
    out = ROOT / "search13_eigenrow_givens_seven_fit_result.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": out.name,
                      "best_cycle_error": result["best_cycle_error"]},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
