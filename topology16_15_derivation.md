# Exact 15-CNOT diagonalizer of the four-qubit right shift

The chronological gate list is [`topology16_15_exact_candidate.json`](topology16_15_exact_candidate.json). Qubit 0 is the most significant bit, and $V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle$. Every one-qubit angle in this gate list is an exact rational multiple of $\pi$. The circuit uses no ancillas.

## Six-CNOT parity prefix

Represent a computational-basis parity by its four-bit mask. The original diagonal parity correction together with the first Fourier CZ has nonlocal (R_Z) parity angles

| Mask | 1011 | 1101 | 1110 | 0110 |
|---|---:|---:|---:|---:|
| Original angle | π/8 | π/4 | 3π/8 | −π/2 |

The four triple masks (1011,1101,1110,0111) form one orbit under the right shift. Multiplying the correction by (R_Z(-\pi/4)) on **each** of these four parity masks commutes with (V_4). The required nonlocal angles become −π/8 on 1011, zero on 1101, π/8 on 1110, −π/4 on 0111, and −π/2 on 0110. The single-wire angles remain (5π/8,3π/4,3π/8) on wires 1,2,3.

The first 13 gates of the candidate realize those angles with six CNOTs. Propagating GF(2) parity masks through the CNOTs gives the final wire masks ((1000,0101,1010,0001)). The chronological inverse CNOTs ((3\to1),(0\to2)) restore the identity mask tuple. Thus this prefix followed by those two inverse CNOTs is the gauged diagonal correction. The subsequent Fourier tail is exactly `baseline_18.json` gates 18 onward; the earlier Fourier CZ has already been absorbed into the parity correction.

## Two exact butterfly replacements

Let (F_{02}) and (F_{13}) be `baseline_18.json` gates 18–27 and 28–37, respectively. They act on disjoint wire pairs, so the inverse CNOT (3\to1) commutes through (F_{02}). The source circuit can therefore be grouped chronologically as

\[
P_6;\ [\operatorname{CX}_{0\to2};F_{02}];\ [\operatorname{CX}_{3\to1};F_{13}];\ T,
\]

where (P_6) is the six-CNOT prefix and (T) is baseline gates 38 onward.

`topology16_15_exact.py` replaces each bracket by a **two-CNOT** circuit. The source and replacement matrices differ only by global phase:

\[
B_{02}=e^{3\pi i/4}F_{02}\operatorname{CX}_{0\to2},\qquad
B_{13}=e^{-\pi i/4}F_{13}\operatorname{CX}_{3\to1}.
\]

Both source blocks have Weyl coordinates ((\pi/4,\pi/8,0)). Their common entangler has the exact two-CNOT identity

\[
R_{XX}(-\pi/2)R_{YY}(-\pi/4)
=K\operatorname{CX}_{c\to t}
 [R_{X,c}(-\pi/2)R_{Z,t}(-\pi/4)]
 \operatorname{CX}_{c\to t}K^\dagger,
\qquad K=R_{X,c}(\pi/2)R_{X,t}(\pi/2).
\]

The one-qubit Weyl factors used around this entangler are explicitly listed in `topology16_exact_block.py` (`exact_block`) and `topology16_15_exact.py` (`reverse_block`). They contain only rational-π rotations. `python3 topology16_15_block_identity.py` proves both 4×4 block identities over (\mathbb Q(\zeta_{32})). The product of their global phases is (i), which does not affect diagonalization.

The count is (6+2+2+5=15) CNOTs: six in (P_6), two in each replacement, and five in (T).

## Exact and independent checks

- `python3 exact_check.py topology16_15_exact_candidate.json` proves (U V_4=D U) over (\mathbb Q(\zeta_{32})), with eigenvalue multiplicities (6,3,4,3) for (1,i,-1,-i).
- Converting every (R_Y(\theta)) to the chronological sequence (R_Z(-\pi/2),R_X(\theta),R_Z(\pi/2)) gives `topology16_15_rxrz_candidate.json`. `python3 delete_search_exact_audit.py topology16_15_rxrz_candidate.json` independently proves the same identity in the integer cyclotomic ring (\mathbb Z[z]/(z^{16}+1)), with a global dyadic denominator.
- Qiskit `Operator` simulation of the original exact-angle candidate gives maximum off-diagonal entry (7.51\times10^{-16}); the local checker gives (5.99\times10^{-16}). Qiskit transpilation retains 15 CNOTs.

The circuit diagonalizes $V_4$. Its transformed total-spin Casimir is not diagonal, so this certificate does not establish a strong Schur transform.
