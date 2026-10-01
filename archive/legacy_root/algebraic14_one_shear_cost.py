"""Compile topology run124's one-shear parity coordinates exactly/RCCX."""

import json
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

from algebraic14_orbit_relative_shears import emit, phase_defects
from algebraic14_orbit_shear_cx import linear_bfs, shortest_form


ROOT = Path(__file__).parent
FORM = (15, 3, 5, 0, 0)


def main():
    action = json.loads((ROOT / "topology14_one_shear_action_result.json").read_text())
    expected = action["coordinate_map"]
    order, parent = linear_bfs()
    conjugator = shortest_form(FORM, order, parent)
    records = []
    for relative in (False, True):
        circuit = QuantumCircuit(4)
        emit(circuit, FORM, conjugator, relative)
        circuit.cx(3, 1)  # physical 0 -> 2
        circuit.cx(2, 1)  # physical 1 -> 2
        phases = phase_defects(circuit, expected)
        compiled = transpile(circuit, basis_gates=["u", "cx"],
                             optimization_level=3, seed_transpiler=42)
        assert np.max(np.abs(Operator(compiled).data -
                             Operator(circuit).data)) < 1e-8
        records.append({
            "relative": relative,
            "linear_conjugator_cx": len(conjugator["linear_gates"]),
            "naive_cx": 2 * len(conjugator["linear_gates"]) +
                        (3 if relative else 6) + 2,
            "compiled_cx": int(compiled.count_ops().get("cx", 0)),
            "phase_exponents_mod8": phases,
            "phase_is_trivial": len(set(phases)) == 1,
        })
    result = {"form": FORM, "expected_coordinate_map": expected,
              "records": records}
    (ROOT / "algebraic14_one_shear_cost_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
