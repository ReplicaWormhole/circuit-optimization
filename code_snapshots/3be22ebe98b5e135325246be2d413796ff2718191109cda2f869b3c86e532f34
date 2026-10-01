"""Search six-CNOT parity prefixes with short nonidentity endpoints.

The diagonal target is the B4-plus-first-CZ parity polynomial from run 6.
Enumerate all reachable six-CNOT prefixes that visit its four nonlocal masks.
For endpoints within three CNOTs of the identity, prepend the inverse linear
map to the nine-CNOT Fourier tail and ask Qiskit to synthesize the combined
tail. Emit the best complete candidate, if any. This is a bounded ansatz only.
"""

import json
from collections import deque
from pathlib import Path

from qiskit import transpile

from algebraic_parity_bfs import ANGLE, SINGLE_ANGLES
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent
BASIS = (8, 4, 2, 1)
REQUIRED = tuple(ANGLE)
FLAGS = {mask: 1 << i for i, mask in enumerate(REQUIRED)}
ALL_FLAGS = (1 << len(REQUIRED)) - 1


def endpoints():
    start = (BASIS, 0)
    queue = deque([(start, 0)])
    previous = {start: None}
    depth_of = {start: 0}
    while queue:
        (rows, flags), depth = queue.popleft()
        if depth == 6:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                nxt = list(rows)
                nxt[target] ^= rows[control]
                nxt = tuple(nxt)
                state = (nxt, flags | FLAGS.get(nxt[target], 0))
                if state not in previous:
                    previous[state] = ((rows, flags), (control, target))
                    depth_of[state] = depth + 1
                    queue.append((state, depth + 1))
    found = []
    for state, depth in depth_of.items():
        if state[1] != ALL_FLAGS or depth != 6:
            continue
        edges = []
        cursor = state
        while cursor != start:
            parent, edge = previous[cursor]
            edges.append(edge)
            cursor = parent
        found.append((state[0], edges[::-1]))
    return found, len(previous)


def inverse_paths(max_depth=3):
    queue = deque([(BASIS, ())])
    found = {BASIS: ()}
    while queue:
        rows, path = queue.popleft()
        if len(path) == max_depth:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                nxt = list(rows)
                nxt[target] ^= rows[control]
                nxt = tuple(nxt)
                if nxt not in found:
                    new_path = path + ((control, target),)
                    found[nxt] = tuple(reversed(new_path))
                    queue.append((nxt, new_path))
    return found


def prefix_gates(edges):
    gates = [{"gate": "rz", "qubit": q, "theta": theta}
             for q, theta in SINGLE_ANGLES.items()]
    rows = list(BASIS)
    seen = set()
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        mask = rows[target]
        if mask in ANGLE and mask not in seen:
            gates.append({"gate": "rz", "qubit": target,
                          "theta": ANGLE[mask]})
            seen.add(mask)
    assert seen == set(REQUIRED)
    return gates


def main():
    ends, states = endpoints()
    inverses = inverse_paths(3)
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    evaluated = 0
    distance_histogram = {}
    best = None
    best_candidate = None
    for rows, edges in ends:
        inverse = inverses.get(rows)
        distance = len(inverse) if inverse is not None else ">3"
        distance_histogram[str(distance)] = distance_histogram.get(str(distance), 0) + 1
        if inverse is None:
            continue
        tail_with_inverse = {"n": 4, "gates": [
            *({"gate": "cx", "control": c, "target": t} for c, t in inverse),
            *tail]}
        optimized = transpile(to_qiskit(tail_with_inverse),
                              basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        candidate = {"n": 4, "gates": prefix_gates(edges)
                     + from_qiskit(optimized)["gates"]}
        check = evaluate(candidate)
        evaluated += 1
        record = {"endpoint_rows": rows, "inverse_cnot_count": len(inverse),
                  "tail_cnot_count": optimized.count_ops().get("cx", 0),
                  "total_cnot_count": check["cnot_count"],
                  "off_diagonal_error": check["off_diagonal_error"],
                  "valid": check["valid_diagonalizer"]}
        if best is None or (record["total_cnot_count"],
                            record["off_diagonal_error"]) < (
                                best["total_cnot_count"], best["off_diagonal_error"]):
            best, best_candidate = record, candidate
    if best_candidate is not None:
        path = ROOT / "algebraic16_to15_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        best["candidate"] = path.name
    print(json.dumps({"states_explored": states, "six_cnot_endpoints": len(ends),
                      "distance_histogram": distance_histogram,
                      "tails_evaluated": evaluated, "best": best}, indent=2))


if __name__ == "__main__":
    main()
