# Candidate gate lists

Circuit JSON lists gates in chronological order. Qubit zero is the most
significant bit; there are no implicit wire swaps. For example:

```json
{
  "n": 4,
  "gates": [
    {"gate": "h", "qubit": 0},
    {"gate": "cx", "control": 0, "target": 3},
    {"gate": "u3", "qubit": 1, "theta": "pi/2", "phi": 0, "lam": "-pi/4"}
  ]
}
```

The one-qubit gates are `h`, `x`, `s`, `sdg`, `rx`, `ry`, `rz`, and `u3`.
Angles may be real numbers in radians or simple multiples of pi, such as
`-3*pi/4`. The two-qubit gate in this format is `cx`. The example illustrates
the format; it is not a diagonalizer.

From the repository root, run `python3 check_circuit.py candidate.json`.
The checker reports the CNOT count, unitarity error, maximum off-diagonal
entry of `U V4 U†`, fourth-root eigenvalue error, output labels, and a
total-spin diagnostic. It permits arbitrary eigenvalue order and orthonormal
bases within degenerate eigenspaces. A passing numerical check does not prove
exactness. Exact claims require a complete gate specification and an exact
certificate; see the [incumbent index](../collaboration/incumbents.json).
