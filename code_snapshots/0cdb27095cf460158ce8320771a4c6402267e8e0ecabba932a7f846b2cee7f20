"""Estimate CNOT cost of exact two-shear orbit preprocessors.

For each quadratic shear, find a shortest 4-bit CNOT network L satisfying
S=L^-1 CCX L with possibly complemented controls. The CNOT network search
is exhaustive over GL(4,2). For every parallel-plane two-shear solution,
rank by the resulting unoptimized cost; Qiskit then compiles the eight
best exact permutation circuits to U/CX as a heuristic joint estimate.
"""

import json
from collections import deque
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

from algebraic14_orbit_shear import orbits, parity, plane_direction
from algebraic14_orbit_two_shears import unique_shears


ROOT = Path(__file__).parent
EDGES = [(c, t) for c in range(4) for t in range(4) if c != t]
IDENTITY = (8, 4, 2, 1)


def update(rows, control, target):
    result = list(rows)
    result[target] ^= result[control]
    return tuple(result)


def linear_bfs():
    parent = {IDENTITY: None}
    queue = deque([IDENTITY])
    order = []
    while queue:
        rows = queue.popleft()
        order.append(rows)
        for edge in EDGES:
            child = update(rows, *edge)
            if child in parent:
                continue
            parent[child] = (rows, edge)
            queue.append(child)
    assert len(order) == 20160
    return order, parent


def shortest_form(form, order, parent):
    v, l1, l2, _, _ = form
    for rows in order:
        if l1 not in rows or l2 not in rows:
            continue
        target_rows = [j for j, mask in enumerate(rows) if parity(mask & v)]
        if len(target_rows) != 1:
            continue
        target = target_rows[0]
        c1, c2 = rows.index(l1), rows.index(l2)
        if target in (c1, c2):
            continue
        gates = []
        state = rows
        while parent[state] is not None:
            previous, edge = parent[state]
            gates.append(edge)
            state = previous
        gates.reverse()
        return {"rows": rows, "control_wires": (c1, c2),
                "target_wire": target, "linear_gates": gates}
    raise RuntimeError(f"No linear conjugator for {form}")


def emit_shear(circuit, form, representation):
    _, _, _, b1, b2 = form
    c1, c2 = representation["control_wires"]
    target = representation["target_wire"]
    gates = representation["linear_gates"]
    for c, t in gates:
        circuit.cx(3 - c, 3 - t)
    if b1:
        circuit.x(3 - c1)
    if b2:
        circuit.x(3 - c2)
    circuit.ccx(3 - c1, 3 - c2, 3 - target)
    if b1:
        circuit.x(3 - c1)
    if b2:
        circuit.x(3 - c2)
    for c, t in reversed(gates):
        circuit.cx(3 - c, 3 - t)


def permutation(circuit):
    matrix = Operator(circuit).data
    output = []
    for column in range(16):
        row = int(np.argmax(abs(matrix[:, column])))
        assert abs(matrix[row, column] - 1) < 1e-9
        assert np.count_nonzero(abs(matrix[:, column]) > 1e-9) == 1
        output.append(row)
    return output


def main():
    order, parent = linear_bfs()
    representatives = unique_shears()
    items = list(representatives.items())
    cycles = [cycle for cycle in orbits() if len(cycle) == 4]
    solutions = []
    for first_perm, first_form in items:
        intermediate = [[first_perm[x] for x in cycle] for cycle in cycles]
        for second_perm, second_form in items:
            transformed = [[second_perm[x] for x in cycle]
                           for cycle in intermediate]
            directions = [plane_direction(cycle)
                          for cycle in transformed]
            if directions[0] is None or not (
                    directions[0] == directions[1] == directions[2]):
                continue
            solutions.append((first_perm, first_form,
                              second_perm, second_form))
    assert len(solutions) == 128
    forms = {form for _, first, _, second in solutions
             for form in (first, second)}
    shortest = {form: shortest_form(form, order, parent) for form in forms}
    records = []
    for first_perm, first_form, second_perm, second_form in solutions:
        first = shortest[first_form]
        second = shortest[second_form]
        raw_cx = 2 * len(first["linear_gates"]) + 2 * len(
            second["linear_gates"]) + 12
        records.append({"first_form": first_form,
                        "second_form": second_form,
                        "first_perm": first_perm,
                        "second_perm": second_perm,
                        "first_linear_cx": len(first["linear_gates"]),
                        "second_linear_cx": len(second["linear_gates"]),
                        "raw_cx_estimate": raw_cx})
    records.sort(key=lambda record: (record["raw_cx_estimate"],
                                     record["first_form"],
                                     record["second_form"]))
    checked = []
    for record in records[:8]:
        circuit = QuantumCircuit(4)
        emit_shear(circuit, record["first_form"],
                   shortest[tuple(record["first_form"])])
        emit_shear(circuit, record["second_form"],
                   shortest[tuple(record["second_form"])])
        expected = [record["second_perm"][record["first_perm"][x]]
                    for x in range(16)]
        assert permutation(circuit) == expected
        compiled = transpile(circuit, basis_gates=["u", "cx"],
                             optimization_level=3, seed_transpiler=42)
        assert np.max(np.abs(Operator(compiled).data -
                             Operator(circuit).data)) < 1e-8
        checked.append({"first_form": record["first_form"],
                        "second_form": record["second_form"],
                        "raw_cx_estimate": record["raw_cx_estimate"],
                        "compiled_cx": int(compiled.count_ops().get("cx", 0)),
                        "permutation": expected})
    result = {"gl4_states": len(order),
              "parallel_plane_solutions": len(records),
              "distinct_shears_used": len(forms),
              "min_raw_cx_estimate": records[0]["raw_cx_estimate"],
              "top_raw_records": [{key: record[key] for key in
                                   ("first_form", "second_form",
                                    "first_linear_cx", "second_linear_cx",
                                    "raw_cx_estimate")}
                                  for record in records[:8]],
              "compiled_top_eight": checked}
    (ROOT / "algebraic14_orbit_shear_cx_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("gl4_states", "parallel_plane_solutions",
                       "distinct_shears_used", "min_raw_cx_estimate",
                       "top_raw_records", "compiled_top_eight")}, indent=2))


if __name__ == "__main__":
    main()
