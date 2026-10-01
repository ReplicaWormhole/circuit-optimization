"""Enumerate rotation-invariant triple-parity gauges at the 17-CNOT boundary.

The four three-body Z masks form one orbit under the right cyclic shift.
Adding the same parity rotation to all four commutes with the shift. Choose
the four offsets that cancel one of the original three-body rotations, then
find the shortest closed CNOT parity network for the remaining three masks
and the pair mask from the first Fourier CZ.
"""

import json
from pathlib import Path

from algebraic_parity_bfs import solve
from check_circuit import evaluate


ROOT = Path(__file__).parent
TRIPLES = (0b1011, 0b1101, 0b1110, 0b0111)
COEFFICIENTS = (1, 2, 3, 0)  # multiples of pi/8
PAIR = 0b0110
SINGLES = {1: "5*pi/8", 2: "3*pi/4", 3: "3*pi/8"}


def angle(coefficient):
    if coefficient == 0:
        raise ValueError("zero phase")
    return f"{coefficient}*pi/8"


def emit(edges, phases):
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    gates = [{"gate": "rz", "qubit": q, "theta": t}
             for q, t in SINGLES.items()]
    rows = [8, 4, 2, 1]
    seen = set()
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        mask = rows[target]
        if mask in phases and mask not in seen:
            gates.append({"gate": "rz", "qubit": target,
                          "theta": phases[mask]})
            seen.add(mask)
    assert rows == [8, 4, 2, 1]
    assert seen == set(phases)
    # Skip B4 (first 15 gates) and the first Fourier CZ (next 3 gates).
    assert baseline["gates"][15:18] == [
        {"gate": "h", "qubit": 1},
        {"gate": "cx", "control": 2, "target": 1},
        {"gate": "h", "qubit": 1},
    ]
    gates.extend(baseline["gates"][18:])
    return {"n": 4, "gates": gates}


def main():
    results = []
    best = None
    for cancelled in range(4):
        shift = -COEFFICIENTS[cancelled]
        phases = {mask: angle(coefficient + shift)
                  for mask, coefficient in zip(TRIPLES, COEFFICIENTS)
                  if coefficient + shift}
        phases[PAIR] = "-pi/2"
        edges, states, transitions = solve(tuple(phases))
        candidate = emit(edges, phases)
        evaluation = evaluate(candidate)
        path = ROOT / f"algebraic16_gauge_{cancelled}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        record = {"cancelled_triple_mask": TRIPLES[cancelled],
                  "shift_pi_over_8": shift,
                  "required_masks": list(phases),
                  "prefix_cnot_count": len(edges),
                  "total_cnot_count": evaluation["cnot_count"],
                  "off_diagonal_error": evaluation["off_diagonal_error"],
                  "states_explored": states,
                  "transitions": transitions,
                  "candidate": path.name}
        results.append(record)
        if best is None or record["total_cnot_count"] < best["total_cnot_count"]:
            best = record
    print(json.dumps({"results": results, "best": best}, indent=2))


if __name__ == "__main__":
    main()
