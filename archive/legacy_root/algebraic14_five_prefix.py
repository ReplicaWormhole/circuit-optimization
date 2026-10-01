"""Bounded search for five-CNOT parity prefixes and absorbable endpoints.

The three-body parity masks form one cyclic orbit. For each of the four
invariant gauge offsets that removes one mask, enumerate five-CNOT prefixes
visiting all required nonlocal masks. For endpoints within four CNOTs of the
identity, synthesize the inverse endpoint with the Fourier tail using Qiskit
level 3, seed 42. This only explores this parity-network ansatz.
"""

import json
from collections import deque
from pathlib import Path

from qiskit import transpile

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent
BASIS = (8, 4, 2, 1)
TRIPLES = (11, 13, 14, 7)
COEFFICIENTS = (1, 2, 3, 0)
SINGLES = {1: "5*pi/8", 2: "3*pi/4", 3: "3*pi/8"}
PAIR = 6


def phase_map(cancelled):
    shift = -COEFFICIENTS[cancelled]
    phases = {mask: f"{coefficient + shift}*pi/8"
              for mask, coefficient in zip(TRIPLES, COEFFICIENTS)
              if coefficient + shift}
    phases[PAIR] = "-pi/2"
    return phases


def enumerate_prefixes(required, depth_limit=5):
    flags = {mask: 1 << i for i, mask in enumerate(required)}
    all_flags = (1 << len(required)) - 1
    start = (BASIS, 0)
    queue = deque([(start, 0)])
    previous = {start: None}
    distances = {start: 0}
    while queue:
        (rows, seen), depth = queue.popleft()
        if depth == depth_limit:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                new_rows = list(rows)
                new_rows[target] ^= rows[control]
                new_rows = tuple(new_rows)
                state = (new_rows, seen | flags.get(new_rows[target], 0))
                if state not in previous:
                    previous[state] = ((rows, seen), (control, target))
                    distances[state] = depth + 1
                    queue.append((state, depth + 1))
    endpoints = []
    for state, depth in distances.items():
        if state[1] != all_flags or depth != depth_limit:
            continue
        edges = []
        cursor = state
        while cursor != start:
            parent, edge = previous[cursor]
            edges.append(edge)
            cursor = parent
        endpoints.append((state[0], edges[::-1]))
    return endpoints, len(previous)


def inverse_paths(max_depth=4):
    queue = deque([(BASIS, ())])
    result = {BASIS: ()}
    while queue:
        rows, edges = queue.popleft()
        if len(edges) == max_depth:
            continue
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                nxt = list(rows)
                nxt[target] ^= rows[control]
                nxt = tuple(nxt)
                if nxt not in result:
                    path = edges + ((control, target),)
                    result[nxt] = tuple(reversed(path))
                    queue.append((nxt, path))
    return result


def prefix_gates(edges, phases):
    gates = [{"gate": "rz", "qubit": q, "theta": theta}
             for q, theta in SINGLES.items()]
    rows = list(BASIS)
    seen = set()
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        mask = rows[target]
        if mask in phases and mask not in seen:
            gates.append({"gate": "rz", "qubit": target,
                          "theta": phases[mask]})
            seen.add(mask)
    assert seen == set(phases)
    return gates


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    inverses = inverse_paths(4)
    results = []
    overall = None
    overall_candidate = None
    for cancelled in range(4):
        phases = phase_map(cancelled)
        ends, states = enumerate_prefixes(tuple(phases))
        histogram = {}
        evaluated = 0
        best = None
        best_candidate = None
        for rows, edges in ends:
            inverse = inverses.get(rows)
            distance = len(inverse) if inverse is not None else ">4"
            histogram[str(distance)] = histogram.get(str(distance), 0) + 1
            if inverse is None:
                continue
            modified_tail = {"n": 4, "gates": [
                *({"gate": "cx", "control": c, "target": t}
                  for c, t in inverse), *tail]}
            optimized = transpile(to_qiskit(modified_tail),
                                  basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = {"n": 4, "gates": prefix_gates(edges, phases)
                         + from_qiskit(optimized)["gates"]}
            check = evaluate(candidate)
            evaluated += 1
            record = {"endpoint_rows": rows,
                      "inverse_distance": len(inverse),
                      "tail_cnot_count": optimized.count_ops().get("cx", 0),
                      "total_cnot_count": check["cnot_count"],
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid": check["valid_diagonalizer"]}
            if best is None or (record["total_cnot_count"],
                                record["off_diagonal_error"]) < (
                                    best["total_cnot_count"], best["off_diagonal_error"]):
                best, best_candidate = record, candidate
        if best_candidate is not None:
            path = ROOT / f"algebraic14_gauge_{cancelled}_best.json"
            path.write_text(json.dumps(best_candidate, indent=2) + "\n")
            best["candidate"] = path.name
            if overall is None or (best["total_cnot_count"],
                                   best["off_diagonal_error"]) < (
                                       overall["total_cnot_count"],
                                       overall["off_diagonal_error"]):
                overall, overall_candidate = best, path.name
        results.append({"cancelled_triple_mask": TRIPLES[cancelled],
                        "required_masks": list(phases),
                        "states_explored": states,
                        "five_cnot_endpoints": len(ends),
                        "inverse_distance_histogram": histogram,
                        "tails_evaluated": evaluated,
                        "best": best})
    print(json.dumps({"results": results,
                      "overall_best": overall,
                      "overall_candidate": overall_candidate}, indent=2))


if __name__ == "__main__":
    main()
