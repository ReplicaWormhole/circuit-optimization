"""Search seven-CNOT parity prefixes ending one CNOT away from identity.

For an endpoint C and diagonal target D, the prefix implements C D. The
Fourier tail must then implement F_tail C. Qiskit is used only to test whether
that modified tail can be transpiled with its original nine CNOTs. Every
complete gate list, if one is found, is checked against V4 afterwards.
"""

import json
from collections import deque
from pathlib import Path

from qiskit import transpile

from algebraic_parity_bfs import ANGLE, SINGLE_ANGLES
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent
REQUIRED = tuple(ANGLE)
MASK_FLAGS = {mask: 1 << i for i, mask in enumerate(REQUIRED)}
BASIS = (8, 4, 2, 1)
ALL_FLAGS = (1 << len(REQUIRED)) - 1


def search_prefixes():
    start = (BASIS, 0)
    queue = deque([(start, 0)])
    previous = {start: None}
    depth_of = {start: 0}
    while queue:
        (rows, flags), depth = queue.popleft()
        if depth == 7:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                new_rows = list(rows)
                new_rows[target] ^= rows[control]
                new_rows = tuple(new_rows)
                new_flags = flags | MASK_FLAGS.get(new_rows[target], 0)
                state = (new_rows, new_flags)
                if state not in previous:
                    previous[state] = ((rows, flags), (control, target))
                    depth_of[state] = depth + 1
                    queue.append((state, depth + 1))
    found = []
    for control in range(4):
        for target in range(4):
            if control == target:
                continue
            rows = list(BASIS)
            rows[target] ^= rows[control]
            state = (tuple(rows), ALL_FLAGS)
            if state not in previous:
                continue
            edges = []
            current = state
            while current != start:
                parent, edge = previous[current]
                edges.append(edge)
                current = parent
            found.append((control, target, edges[::-1]))
    return found, len(previous)


def prefix_gates(edges):
    gates = [{"gate": "rz", "qubit": q, "theta": theta}
             for q, theta in SINGLE_ANGLES.items()]
    rows = list(BASIS)
    seen = set()
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        if rows[target] in ANGLE and rows[target] not in seen:
            gates.append({"gate": "rz", "qubit": target,
                          "theta": ANGLE[rows[target]]})
            seen.add(rows[target])
    assert seen == set(REQUIRED)
    return gates


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    prefixes, states = search_prefixes()
    results = []
    best = None
    for control, target, edges in prefixes:
        tail_with_c = {"n": 4, "gates": [
            {"gate": "cx", "control": control, "target": target}, *tail]}
        optimized = transpile(to_qiskit(tail_with_c), basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=42)
        tail_candidate = from_qiskit(optimized)
        candidate = {"n": 4, "gates": prefix_gates(edges)
                     + tail_candidate["gates"]}
        check = evaluate(candidate)
        path = ROOT / f"algebraic16_endpoint_{control}_{target}.json"
        path.write_text(json.dumps(candidate, indent=2) + "\n")
        record = {"endpoint_cnot": [control, target],
                  "prefix_cnot_count": len(edges),
                  "tail_cnot_count": optimized.count_ops().get("cx", 0),
                  "total_cnot_count": check["cnot_count"],
                  "off_diagonal_error": check["off_diagonal_error"],
                  "valid": check["valid_diagonalizer"],
                  "candidate": path.name}
        results.append(record)
        if best is None or (record["total_cnot_count"], record["off_diagonal_error"]) < (
                best["total_cnot_count"], best["off_diagonal_error"]):
            best = record
    print(json.dumps({"states": states, "prefixes": len(prefixes),
                      "results": results, "best": best}, indent=2))


if __name__ == "__main__":
    main()
