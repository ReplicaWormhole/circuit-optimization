"""Screen relative-phase Toffoli variants of two-shear orbit preprocessors.

For all 128 parallel-plane two-shear maps, replace either or both exact CCX
gates by three-CX RCCX. A replacement has the correct basis permutation but
introduces input phases. Those phases must be constant along each original
V4 orbit to preserve the same transformed cycle action for a shared Fourier
block. Report the phase obstruction and Qiskit CX counts of selected cases.
"""

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

from algebraic14_orbit_shear import orbits, plane_direction
from algebraic14_orbit_shear_cx import linear_bfs, shortest_form
from algebraic14_orbit_two_shears import unique_shears


ROOT = Path(__file__).parent
MODES = ((True, True), (True, False), (False, True))


def emit(circuit, form, representation, relative):
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
    if relative:
        circuit.rccx(3 - c1, 3 - c2, 3 - target)
    else:
        circuit.ccx(3 - c1, 3 - c2, 3 - target)
    if b1:
        circuit.x(3 - c1)
    if b2:
        circuit.x(3 - c2)
    for c, t in reversed(gates):
        circuit.cx(3 - c, 3 - t)


def phase_defects(circuit, expected):
    matrix = Operator(circuit).data
    phases = []
    for x in range(16):
        output = expected[x]
        phase = matrix[output, x]
        assert abs(abs(phase) - 1) < 1e-9
        assert np.count_nonzero(abs(matrix[:, x]) > 1e-9) == 1
        angle = math.atan2(phase.imag, phase.real)
        k = round(4 * angle / math.pi) % 8
        assert abs(phase - np.exp(1j * k * math.pi / 4)) < 1e-9
        phases.append(k)
    return phases


def main():
    order, parent = linear_bfs()
    items = list(unique_shears().items())
    cycles = [cycle for cycle in orbits() if len(cycle) == 4]
    all_cycles = orbits()
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
        expected = [second_perm[first_perm[x]] for x in range(16)]
        for mode in MODES:
            circuit = QuantumCircuit(4)
            emit(circuit, first_form, shortest[first_form], mode[0])
            emit(circuit, second_form, shortest[second_form], mode[1])
            phases = phase_defects(circuit, expected)
            orbit_constant = all(
                len({phases[x] for x in cycle}) == 1
                for cycle in all_cycles)
            varying_cycles = [cycle for cycle in all_cycles
                              if len({phases[x] for x in cycle}) > 1]
            first_depth = len(shortest[first_form]["linear_gates"])
            second_depth = len(shortest[second_form]["linear_gates"])
            raw_cx = (2 * first_depth + 2 * second_depth +
                      sum(3 if bit else 6 for bit in mode))
            records.append({"first_form": first_form,
                            "second_form": second_form,
                            "mode": mode,
                            "raw_cx_estimate": raw_cx,
                            "phase_exponents_mod8": phases,
                            "orbit_constant": orbit_constant,
                            "varying_orbit_count": len(varying_cycles)})
    records.sort(key=lambda record: (
        0 if record["orbit_constant"] else 1,
        record["raw_cx_estimate"],
        record["varying_orbit_count"],
        record["first_form"], record["second_form"], record["mode"]))
    compile_records = records[:12]
    for record in compile_records:
        circuit = QuantumCircuit(4)
        first_form = tuple(record["first_form"])
        second_form = tuple(record["second_form"])
        emit(circuit, first_form, shortest[first_form], record["mode"][0])
        emit(circuit, second_form, shortest[second_form], record["mode"][1])
        optimized = transpile(circuit, basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        assert np.max(np.abs(Operator(optimized).data -
                             Operator(circuit).data)) < 1e-8
        record["compiled_cx"] = int(optimized.count_ops().get("cx", 0))
    by_mode = {}
    for mode in MODES:
        subset = [record for record in records if record["mode"] == mode]
        by_mode["".join("R" if bit else "E" for bit in mode)] = {
            "trials": len(subset),
            "orbit_constant": sum(record["orbit_constant"]
                                  for record in subset),
            "minimum_raw_cx": min(record["raw_cx_estimate"]
                                  for record in subset),
            "minimum_varying_orbits": min(record["varying_orbit_count"]
                                         for record in subset)}
    result = {"total_trials": len(records),
              "phase_exponents_in_units_pi_over4": True,
              "by_mode": by_mode,
              "common_phase_patterns": [
                  {"exponents_mod8": list(pattern), "count": count}
                  for pattern, count in Counter(
                      tuple(record["phase_exponents_mod8"])
                      for record in records).most_common(8)],
              "compiled_best_twelve": compile_records}
    (ROOT / "algebraic14_orbit_relative_shears_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("total_trials", "by_mode", "compiled_best_twelve")},
                     indent=2))


if __name__ == "__main__":
    main()
