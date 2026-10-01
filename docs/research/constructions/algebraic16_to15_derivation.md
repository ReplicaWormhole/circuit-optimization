# Exact 15-CNOT four-cycle diagonalizer

The chronological gate list is `topology16_15_exact_candidate.json`; all its
one-qubit angles are rational multiples of \(\pi\). Qubit 0 is the most
significant bit. This diagonalizes the right cyclic shift \(V_4\), with no
ancilla and 15 CNOTs. It does not diagonalize total spin.

Let \(B_4\) be the paper's diagonal parity correction and
\(P=Z_0Z_1Z_2Z_3\). Remove the full-parity rotation
\(R_P(-3\pi/4)\) from \(B_4\), obtaining
\(B'_4=B_4 R_P(3\pi/4)\). Since \([P,V_4]=0\), this preserves
diagonalization. Next define

\[
G=\prod_{M\in\{Z_0Z_2Z_3,\,Z_0Z_1Z_3,\,Z_0Z_1Z_2,\,Z_1Z_2Z_3\}}
R_M(-\pi/4).
\]

The shift permutes these four three-body masks, so \([G,V_4]=0\). Thus
\(F_4 B'_4 G\) diagonalizes \(V_4\) exactly. The first Fourier CZ is combined
with the diagonal correction. Its resulting nonlocal parity rotations are

| mask on qubits 0–3 | angle |
| --- | ---: |
| 0111 | \(-\pi/4\) |
| 0110 | \(-\pi/2\) |
| 1110 | \(\pi/8\) |
| 1011 | \(-\pi/8\) |

The six-CNOT prefix in the gate list visits these masks in table order, with
local \(R_z\) angles \(5\pi/8,3\pi/4,3\pi/8\) on qubits 1, 2, 3. Its wire
parity rows at the end are `(1000,0101,1010,0001)`, so the prefix implements
\(C\,CZ_{1,2}B'_4G\) up to global phase, where
\(C=CNOT_{3\to1}CNOT_{0\to2}\). The two CNOTs in \(C\) act on disjoint
pairs and commute.

The Fourier tail starts with lifted Hadamards on pairs (0,2) and (1,3).
Compose the first with \(CNOT_{0\to2}\) and the second with
\(CNOT_{3\to1}\). `topology16_15_exact.py` supplies exact two-CNOT
decompositions for both resulting two-qubit unitaries. The remaining Fourier
tail is unchanged. Therefore the entire circuit uses six CNOTs in the open
parity prefix, two in each modified butterfly, and five later: \(6+2+2+5=15\).

`exact_check.py topology16_15_exact_candidate.json` and the independent
`algebraic16_to15_exact_check.check` verify \(U V_4=D U\) over
\(\mathbb Q(\zeta_{32})\) with exact arithmetic. Both return CNOT count 15,
eigenvalue labels `[0,3,1,1,2,2,0,2,0,0,2,0,3,1,3,0]`, and multiplicities
`[6,3,4,3]`. Numerical simulation and Qiskit give residuals near machine
precision. These checks establish an exact diagonalizer; they do not assert
optimality or a stronger Schur transform.
