"""Find short closed CNOT parity networks for B4 composed with the first CZ.

Each wire carries a GF(2)^4 parity mask. A CNOT replaces its target row by
the XOR of target and control. A diagonal Rz on a wire can apply the phase
associated with its current mask. BFS seeks a shortest network that visits
every nonlocal mask needed and returns the row basis to the identity.
"""

import argparse
import json
from collections import deque
from pathlib import Path

from check_circuit import evaluate


REQUIRED = (0b1011, 0b1101, 0b1110, 0b0110)  # 023, 013, 012, 12
ANGLE = {0b1011: "pi/8", 0b1101: "pi/4", 0b1110: "3*pi/8", 0b0110: "-pi/2"}
SINGLE_ANGLES = {1: "5*pi/8", 2: "3*pi/4", 3: "3*pi/8"}


def solve(required=REQUIRED):
    basis = (8, 4, 2, 1)
    visited = tuple(1 << i for i, mask in enumerate(required))
    initial = (basis, 0)
    target = (basis, (1 << len(required)) - 1)
    queue = deque([initial])
    previous = {initial: None}
    transitions = 0
    while queue:
        state = queue.popleft()
        rows, flags = state
        if state == target:
            break
        for control in range(4):
            for target_wire in range(4):
                if control == target_wire:
                    continue
                new_rows = list(rows)
                new_rows[target_wire] ^= rows[control]
                new_rows = tuple(new_rows)
                new_flags = flags
                if new_rows[target_wire] in required:
                    new_flags |= 1 << required.index(new_rows[target_wire])
                new_state = (new_rows, new_flags)
                transitions += 1
                if new_state not in previous:
                    previous[new_state] = (state, (control, target_wire))
                    queue.append(new_state)
    if target not in previous:
        raise RuntimeError("no closed parity network")
    edges = []
    state = target
    while state != initial:
        parent, edge = previous[state]
        edges.append(edge)
        state = parent
    edges.reverse()
    return edges, len(previous), transitions


def emit_candidate(edges, baseline):
    gates = []
    for q, theta in SINGLE_ANGLES.items():
        gates.append({"gate": "rz", "qubit": q, "theta": theta})
    rows = [8, 4, 2, 1]
    seen = set()
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        mask = rows[target]
        if mask in ANGLE and mask not in seen:
            gates.append({"gate": "rz", "qubit": target, "theta": ANGLE[mask]})
            seen.add(mask)
    assert tuple(rows) == (8, 4, 2, 1)
    assert seen == set(ANGLE)
    # The baseline's first Fourier CZ is h(1), cx(2,1), h(1).
    first_cz = next(i for i, gate in enumerate(baseline["gates"])
                    if gate["gate"] == "h" and gate.get("qubit") == 1)
    assert baseline["gates"][first_cz:first_cz+3] == [
        {"gate": "h", "qubit": 1},
        {"gate": "cx", "control": 2, "target": 1},
        {"gate": "h", "qubit": 1},
    ]
    gates.extend(baseline["gates"][first_cz+3:])
    return {"n": 4, "gates": gates}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    edges, states, transitions = solve()
    baseline = json.loads(Path(__file__).with_name("baseline_18.json").read_text())
    candidate = emit_candidate(edges, baseline)
    result = evaluate(candidate)
    if args.output:
        args.output.write_text(json.dumps(candidate, indent=2) + "\n")
    print(json.dumps({"edges": edges, "network_cnot_count": len(edges),
                      "states_explored": states, "transitions": transitions,
                      "candidate_evaluation": result}, indent=2))


if __name__ == "__main__":
    main()
