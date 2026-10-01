"""Definition-level construction of ONE fixed permutation, not a search."""
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
P = [0, 4, 7, 8, 6, 2, 11, 12, 5, 9, 3, 13, 10, 14, 15, 1]

def bit(x, q):
    return (x >> (3-q)) & 1

def anf(q):
    coefficients = [bit(P[x], q) for x in range(16)]
    for b in (1, 2, 4, 8):
        for mask in range(16):
            if mask & b:
                coefficients[mask] ^= coefficients[mask ^ b]
    masks = [m for m, value in enumerate(coefficients) if value]
    for x in range(16):
        assert sum((x & m) == m for m in masks) % 2 == bit(P[x], q)
    return masks

visited = set()
cycles = []
gates = []
transpositions = []
for start in range(16):
    if start in visited:
        continue
    cycle = []
    x = start
    while x not in visited:
        cycle.append(x)
        visited.add(x)
        x = P[x]
    assert x == start
    cycles.append(cycle)
    for end in cycle[1:]:
        path = [start]
        for q in range(4):
            if bit(start, q) != bit(end, q):
                path.append(path[-1] ^ (1 << (3-q)))
        edges = list(zip(path, path[1:]))
        for a, b in edges + edges[-2::-1]:
            target = next(q for q in range(4) if bit(a, q) != bit(b, q))
            gates.append({"gate": "mcx", "target": target,
                          "controls": [{"qubit": q, "value": bit(a, q)}
                                       for q in range(4) if q != target]})
        transpositions.append({"a": start, "b": end, "gray_path": path})

actual = []
for x in range(16):
    y = x
    for g in gates:
        if all(bit(y, c["qubit"]) == c["value"] for c in g["controls"]):
            y ^= 1 << (3-g["target"])
    actual.append(y)
assert actual == P
result = {"fixed_mapping": P, "cycles": cycles,
          "anf_masks_q0_through_q3": [anf(q) for q in range(4)],
          "transpositions": transpositions, "gates": gates,
          "native_three_control_x_count": len(gates),
          "compiled_cx_upper_bound": 34*len(gates),
          "verification": "All16 direct gate evaluations equal fixed mapping; all64 output bits equal ANF",
          "scope": "One deterministic construction, no search or optimized cost claim"}
(OUT / "encoder.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({k: v for k, v in result.items() if k not in ("gates", "transpositions")}))
