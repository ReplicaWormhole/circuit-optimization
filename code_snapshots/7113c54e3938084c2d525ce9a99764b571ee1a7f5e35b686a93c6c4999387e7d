"""Screen simultaneous three-orbit Householder gauges in a pure orbit DFT basis.

Construct exact orbit eigenvectors for V4 (three length-4 orbits, one
length-2 orbit, two fixed states). Match them to exact15 output slots within
each eigenvalue by maximum overlap, then apply rational 3x3 Householder
reflections to the three length-4 rows in selected eigenspaces. Measure
rank signatures after the exact15 six-CX prefix as a filter for <=14
synthesis. This numerical screen does not establish exact rank.
"""

import json
from pathlib import Path

import numpy as np
from qiskit.quantum_info import Operator
from scipy.optimize import linear_sum_assignment

from algebraic14_orbit_shear import orbits, rotate
from algebraic14_threeway_gauge_rank import graph_bound, gauge, ranks
from exact_check import exact_eigenvalue_labels
from qiskit_crosscheck import to_qiskit


ROOT = Path(__file__).parent
ROOTS = (1, 1j, -1, -1j)
VECTORS = ((1, 1, 1), (1, 1, 2), (1, 2, 2))
SCOPES = ((0, 1, 2, 3), (0, 2), (1, 3), (0,), (2,))


def orbit_vectors():
    records = []
    four_cycles = [cycle for cycle in orbits() if len(cycle) == 4]
    for orbit_index, cycle in enumerate(four_cycles):
        for label, root in enumerate(ROOTS):
            row = np.zeros(16, dtype=complex)
            for j, x in enumerate(cycle):
                row[x] = root ** j / 2
            records.append({"label": label, "orbit_index": orbit_index,
                            "kind": "length4", "row": row})
    two = [cycle for cycle in orbits() if len(cycle) == 2]
    assert len(two) == 1
    for label, root in ((0, 1), (2, -1)):
        row = np.zeros(16, dtype=complex)
        for j, x in enumerate(two[0]):
            row[x] = root ** j / np.sqrt(2)
        records.append({"label": label, "orbit_index": None,
                        "kind": "length2", "row": row})
    for x in (0, 15):
        row = np.zeros(16, dtype=complex)
        row[x] = 1
        records.append({"label": 0, "orbit_index": None,
                        "kind": "fixed", "row": row})
    assert len(records) == 16
    return records


def matched_orbit_basis(exact15, exact15_labels):
    records = orbit_vectors()
    u = np.zeros((16, 16), dtype=complex)
    rows_by_label = {label: {} for label in range(4)}
    assignments = []
    for label in range(4):
        output_slots = [j for j, value in enumerate(exact15_labels)
                        if value == label]
        candidates = [record for record in records
                      if record["label"] == label]
        overlap = np.array([
            [abs(np.vdot(exact15[slot], record["row"])) ** 2
             for record in candidates] for slot in output_slots])
        selected_rows, selected_cols = linear_sum_assignment(-overlap)
        for a, b in zip(selected_rows, selected_cols):
            slot = output_slots[a]
            record = candidates[b]
            u[slot] = record["row"]
            assignments.append({"output_row": slot,
                                "label": label,
                                "kind": record["kind"],
                                "orbit_index": record["orbit_index"],
                                "overlap_squared": float(overlap[a, b])})
            if record["kind"] == "length4":
                rows_by_label[label][record["orbit_index"]] = slot
    assert all(len(group) == 3 for group in rows_by_label.values())
    return u, {label: [rows_by_label[label][j] for j in range(3)]
               for label in range(4)}, assignments


def main():
    candidate = json.loads(
        (ROOT / "topology16_15_exact_candidate.json").read_text())
    exact15 = Operator(to_qiskit(candidate)).data
    exact15_labels = exact_eigenvalue_labels(candidate)["output_labels"]
    u, rows_by_label, assignments = matched_orbit_basis(
        exact15, exact15_labels)
    prefix = Operator(to_qiskit(
        {"n": 4, "gates": candidate["gates"][:13]})).data
    shift = np.zeros((16, 16), dtype=complex)
    for x in range(16):
        shift[rotate(x), x] = 1
    assert np.max(np.abs(u @ u.conj().T - np.eye(16))) < 1e-12
    baseline_tail_ranks = ranks(u @ prefix.conj().T)
    results = []
    for vector in VECTORS:
        for scope in SCOPES:
            transformed = gauge(rows_by_label, vector, scope) @ u
            diagonal = transformed @ shift @ transformed.conj().T
            offdiag = float(np.max(np.abs(
                diagonal - np.diag(np.diag(diagonal)))))
            assert offdiag < 1e-12
            tail_ranks = ranks(transformed @ prefix.conj().T)
            results.append({
                "vector": vector, "scope": scope,
                "full_ranks": ranks(transformed),
                "tail_ranks": tail_ranks,
                "tail_cut_graph_bound": graph_bound(tail_ranks),
                "offdiag_error": offdiag,
                "rank_improvement_vs_pure_orbit": sum(
                    tail_ranks[key] < baseline_tail_ranks[key]
                    for key in baseline_tail_ranks),
            })
    results.sort(key=lambda row: (
        row["tail_cut_graph_bound"],
        -row["rank_improvement_vs_pure_orbit"],
        sum(row["tail_ranks"].values()),
        row["vector"], row["scope"]))
    result = {"assignments": assignments,
              "rows_by_label": rows_by_label,
              "pure_orbit_tail_ranks": baseline_tail_ranks,
              "pure_orbit_tail_graph_bound": graph_bound(baseline_tail_ranks),
              "cases": len(results),
              "records": results}
    (ROOT / "algebraic14_pure_orbit_householder_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("rows_by_label", "pure_orbit_tail_ranks",
                       "pure_orbit_tail_graph_bound", "cases")},
                     indent=2))
    print(json.dumps({"top_records": results[:6]}, indent=2))


if __name__ == "__main__":
    main()
