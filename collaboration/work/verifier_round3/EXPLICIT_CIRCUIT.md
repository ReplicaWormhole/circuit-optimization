# Explicit accepted 13-CNOT diagonalizer

Current incumbent, not a proven optimum. The retained bound is 6 <= C_min <= 13.

Gate model: ancilla-free four qubits, arbitrary exact single-qubit rotations and fixed directed CNOT on any ordered pair. Other native entanglers have separate costs.

Qubit 0 is most significant. Execute rows in ascending order; U = G_73 ... G_1. The target is U V4 = D U, where V4|x0 x1 x2 x3> = |x3 x0 x1 x2>. Output eigenvalue order is arbitrary.

R_a(theta) = cos(theta/2) I - i sin(theta/2) sigma_a. The signed half-angle columns specify the exact rotations, including their branches; theta/pi is shown only when supplied exactly. Every radical sqrt denotes the nonnegative real root. Full-angle columns preserve the source metadata. Zero rotations are retained.

Source: `collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`

SHA256: `d500ff2a79c4717ff38d935d71b67bd798b06a7c6899a275fec25570b86b8d13`

Existing acceptance: `collaboration/EXACT13_ACCEPTANCE.json`, certificate runs 370 and 371 using the same full-matrix implementation, with separate relation and scalar audits. This export checks transcription only; it does not repeat or add an exact matrix certificate.

| Step | Gate | Wire(s) | theta/pi | cos(theta/2) | sin(theta/2) | cos(theta) | sin(theta) | Coordinate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | RZ | 0 | 0 | 1 | 0 | 1 | 0 | 0 |
| 2 | RY | 0 | 0 | 1 | 0 | 1 | 0 | 1 |
| 3 | RX | 1 | 3/2 | -sqrt(2)/2 | sqrt(2)/2 | 0 | -1 | 2 |
| 4 | RY | 1 | 0 | 1 | 0 | 1 | 0 | 3 |
| 5 | CX | 0 -> 1 | — | — | — | — | — | — |
| 6 | RX | 2 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 4 |
| 7 | RY | 2 | 1 | 0 | 1 | -1 | 0 | 5 |
| 8 | RZ | 3 | 0 | 1 | 0 | 1 | 0 | 6 |
| 9 | RY | 3 | 0 | 1 | 0 | 1 | 0 | 7 |
| 10 | CX | 3 -> 2 | — | — | — | — | — | — |
| 11 | RZ | 1 | 0 | 1 | 0 | 1 | 0 | 8 |
| 12 | RY | 1 | 0 | 1 | 0 | 1 | 0 | 9 |
| 13 | RX | 2 | 0 | 1 | 0 | 1 | 0 | 10 |
| 14 | RY | 2 | 0 | 1 | 0 | 1 | 0 | 11 |
| 15 | CX | 1 -> 2 | — | — | — | — | — | — |
| 16 | RX | 0 | 7/8 | sqrt(1/2 - sqrt(sqrt(2)/4 + 1/2)/2) | sqrt(sqrt(sqrt(2)/4 + 1/2)/2 + 1/2) | -sqrt(sqrt(2)/4 + 1/2) | sqrt(1/2 - sqrt(2)/4) | 12 |
| 17 | RY | 0 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 13 |
| 18 | RZ | 1 | — | sqrt(1/2 - sqrt(3)/6) | sqrt(sqrt(3)/6 + 1/2) | -sqrt(3)/3 | sqrt(6)/3 | 14 |
| 19 | RY | 1 | — | sqrt(1/2 - sqrt(10)/10) | sqrt(sqrt(10)/10 + 1/2) | -sqrt(10)/5 | sqrt(15)/5 | 15 |
| 20 | CX | 1 -> 0 | — | — | — | — | — | — |
| 21 | RZ | 2 | 0 | 1 | 0 | 1 | 0 | 16 |
| 22 | RY | 2 | 0 | 1 | 0 | 1 | 0 | 17 |
| 23 | RX | 3 | -5/8 | sqrt(1/2 - sqrt(1/2 - sqrt(2)/4)/2) | -sqrt(sqrt(1/2 - sqrt(2)/4)/2 + 1/2) | -sqrt(1/2 - sqrt(2)/4) | -sqrt(sqrt(2)/4 + 1/2) | 18 |
| 24 | RY | 3 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 19 |
| 25 | CX | 2 -> 3 | — | — | — | — | — | — |
| 26 | RX | 1 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 20 |
| 27 | RY | 1 | — | sqrt(1/2 - sqrt(30)/12) | sqrt(sqrt(30)/12 + 1/2) | -sqrt(30)/6 | sqrt(6)/6 | 21 |
| 28 | RZ | 2 | 1/8 | sqrt(sqrt(sqrt(2)/4 + 1/2)/2 + 1/2) | sqrt(1/2 - sqrt(sqrt(2)/4 + 1/2)/2) | sqrt(sqrt(2)/4 + 1/2) | sqrt(1/2 - sqrt(2)/4) | 22 |
| 29 | RY | 2 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 23 |
| 30 | CX | 2 -> 1 | — | — | — | — | — | — |
| 31 | RZ | 0 | 3/2 | -sqrt(2)/2 | sqrt(2)/2 | 0 | -1 | 24 |
| 32 | RY | 0 | 0 | 1 | 0 | 1 | 0 | 25 |
| 33 | RX | 1 | -1/2 | sqrt(2)/2 | -sqrt(2)/2 | 0 | -1 | 26 |
| 34 | RY | 1 | — | sqrt(sqrt(45*sqrt(2)/367 + 115/367)/2 + 1/2) | sqrt(1/2 - sqrt(45*sqrt(2)/367 + 115/367)/2) | sqrt(45*sqrt(2)/367 + 115/367) | sqrt(252/367 - 45*sqrt(2)/367) | 27 |
| 35 | CX | 0 -> 1 | — | — | — | — | — | — |
| 36 | RX | 2 | 0 | 1 | 0 | 1 | 0 | 28 |
| 37 | RY | 2 | 1 | 0 | 1 | -1 | 0 | 29 |
| 38 | RZ | 3 | -1/2 | sqrt(2)/2 | -sqrt(2)/2 | 0 | -1 | 30 |
| 39 | RY | 3 | 0 | 1 | 0 | 1 | 0 | 31 |
| 40 | CX | 3 -> 2 | — | — | — | — | — | — |
| 41 | RZ | 1 | -1/2 | sqrt(2)/2 | -sqrt(2)/2 | 0 | -1 | 32 |
| 42 | RY | 1 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 33 |
| 43 | RX | 2 | 0 | 1 | 0 | 1 | 0 | 34 |
| 44 | RY | 2 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 35 |
| 45 | CX | 1 -> 2 | — | — | — | — | — | — |
| 46 | RX | 0 | 0 | 1 | 0 | 1 | 0 | 36 |
| 47 | RY | 0 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 37 |
| 48 | RZ | 1 | — | sqrt(1/2 - sqrt(sqrt(2)/7 + 4/7)/2) | sqrt(sqrt(sqrt(2)/7 + 4/7)/2 + 1/2) | -sqrt(sqrt(2)/7 + 4/7) | sqrt(3/7 - sqrt(2)/7) | 38 |
| 49 | RY | 1 | — | sqrt(1/2 - sqrt(5*sqrt(2)/48 + 7/12)/2) | sqrt(sqrt(5*sqrt(2)/48 + 7/12)/2 + 1/2) | -sqrt(5*sqrt(2)/48 + 7/12) | sqrt(5/12 - 5*sqrt(2)/48) | 39 |
| 50 | CX | 1 -> 0 | — | — | — | — | — | — |
| 51 | RZ | 2 | -1/2 | sqrt(2)/2 | -sqrt(2)/2 | 0 | -1 | 40 |
| 52 | RY | 2 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 41 |
| 53 | RX | 3 | — | sqrt(1/2 - sqrt(6)/6) | sqrt(sqrt(6)/6 + 1/2) | -sqrt(6)/3 | sqrt(3)/3 | 42 |
| 54 | RY | 3 | 1/3 | sqrt(3)/2 | 1/2 | 1/2 | sqrt(3)/2 | 43 |
| 55 | CX | 2 -> 3 | — | — | — | — | — | — |
| 56 | RX | 1 | — | sqrt(1/2 - sqrt(183/367 - 72*sqrt(2)/367)/2) | sqrt(sqrt(183/367 - 72*sqrt(2)/367)/2 + 1/2) | -sqrt(183/367 - 72*sqrt(2)/367) | sqrt(72*sqrt(2)/367 + 184/367) | 44 |
| 57 | RY | 1 | — | sqrt(sqrt(9/28 - 3*sqrt(2)/28)/2 + 1/2) | sqrt(1/2 - sqrt(9/28 - 3*sqrt(2)/28)/2) | sqrt(9/28 - 3*sqrt(2)/28) | sqrt(3*sqrt(2)/28 + 19/28) | 45 |
| 58 | RZ | 2 | 0 | 1 | 0 | 1 | 0 | 46 |
| 59 | RY | 2 | 1 | 0 | 1 | -1 | 0 | 47 |
| 60 | CX | 2 -> 1 | — | — | — | — | — | — |
| 61 | RZ | 0 | — | -sqrt(1/2 - sqrt(10)/10) | sqrt(sqrt(10)/10 + 1/2) | -sqrt(10)/5 | -sqrt(15)/5 | 48 |
| 62 | RY | 0 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 49 |
| 63 | RX | 1 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 50 |
| 64 | RY | 1 | 1/2 | sqrt(2)/2 | sqrt(2)/2 | 0 | 1 | 51 |
| 65 | CX | 0 -> 1 | — | — | — | — | — | — |
| 66 | RZ | 0 | 0 | 1 | 0 | 1 | 0 | 52 |
| 67 | RY | 0 | 3/4 | sqrt(1/2 - sqrt(2)/4) | sqrt(sqrt(2)/4 + 1/2) | -sqrt(2)/2 | sqrt(2)/2 | 53 |
| 68 | RZ | 1 | — | -sqrt(1/2 - sqrt(6)/6) | -sqrt(sqrt(6)/6 + 1/2) | -sqrt(6)/3 | sqrt(3)/3 | 54 |
| 69 | RY | 1 | 3/4 | sqrt(1/2 - sqrt(2)/4) | sqrt(sqrt(2)/4 + 1/2) | -sqrt(2)/2 | sqrt(2)/2 | 55 |
| 70 | RZ | 2 | 0 | 1 | 0 | 1 | 0 | 56 |
| 71 | RY | 2 | 0 | 1 | 0 | 1 | 0 | 57 |
| 72 | RZ | 3 | — | sqrt(sqrt(15)/10 + 1/2) | -sqrt(1/2 - sqrt(15)/10) | sqrt(15)/5 | -sqrt(10)/5 | 58 |
| 73 | RY | 3 | — | sqrt(1/2 - sqrt(6)/12) | sqrt(sqrt(6)/12 + 1/2) | -sqrt(6)/6 | sqrt(30)/6 | 59 |

D output labels (basis states 0000 through 1111, in binary order). Labels 0, 1, 2, 3 mean 1, i, -1, -i, respectively; D_rr = i**labels[r]:

`[2, 0, 0, 1, 3, 1, 2, 3, 0, 2, 0, 0, 1, 3, 2, 0]`

Reproduce: `python3 collaboration/work/verifier_round3/export_circuit.py`
