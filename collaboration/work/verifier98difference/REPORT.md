# Independent parity-difference Bell-block review

Board 110. This is a by-hand review only. I performed no matrix evaluation,
symbolic computation, search, optimizer run, or fit.

## Conjugation and block action

Let (a=(a_0,a_1)) and (b=(b_0,b_1)) be the Bell labels on
(A=(q_0,q_2)) and (B=(q_1,q_3)), respectively, and put

\[
\sigma_x=(-1)^{x_0x_1},\qquad T|a,b\rangle=\sigma_b|b,a\rangle.
\]

Use the reversible difference relabeling
\[
R|a,b\rangle=|a,d=b\oplus a\rangle.
\]
It is self-inverse. This is exactly the pair of label-bit CNOTs
`CX(q0->q1) CX(q2->q3)`: each target label bit is XORed with the
corresponding control label bit. For an output basis state \(|a,d\rangle\),
\[
R^\dagger|a,d\rangle=|a,a\oplus d\rangle,
\]
so
\[
R T R^\dagger|a,d\rangle
=R\bigl(\sigma_{a\oplus d}|a\oplus d,a\rangle\bigr)
=\sigma_{a\oplus d}|a\oplus d,d\rangle.
\]
Thus d is invariant and the operator on A in sector d is
\[
T_d=CZ_A X_A^d,
\]
where (X_A^d=X_0^{d_0}X_1^{d_1}) acts first. Equivalently, the phase
is applied to the flipped label (a\oplus d). This confirms the stated
orientation and product order.

## Four sector eigenbases

The following pair orderings specify the two-dimensional blocks. In each
ordinary-swap block, (Q=H). In each signed-exchange block with matrix
\(\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}\), the eigenvectors
ordered as (+i,-i) give (Q=SH), so the applied diagonalizer is
\(Q^\dagger=HS^\dagger\), with (S=\operatorname{diag}(1,i)).

| (d) | blocks on A | eigenbasis / eigenvalues |
|---|---|---|
| 00 | singletons 00, 01, 10, 11 | computational basis; (+,+,+,-) |
| 10 | (00,10), (01,11) | first is ordinary swap: (H), eigenvalues (+,-). Second has (T|01\rangle=-|11\rangle, T|11\rangle=|01\rangle): (SH), eigenvalues (+i,-i). |
| 01 | (00,01), (10,11) | first is ordinary swap: (H), eigenvalues (+,-). Second has (T|10\rangle=-|11\rangle, T|11\rangle=|10\rangle): (SH), eigenvalues (+i,-i). |
| 11 | (01,10), (00,11) | first is ordinary swap: (H), eigenvalues (+,-). Second has (T|00\rangle=-|11\rangle, T|11\rangle=|00\rangle): (SH), eigenvalues (+i,-i). |

Together these sectors reproduce multiplicities ((6,4,3,3)) for
((+1,-1,+i,-i)). Since output diagonal order is unrestricted, the
diagonalizer can be written (W=C R), where (C=\bigoplus_d Q_d^\dagger)
is the sector-controlled eigenbasis change. No final (R^\dagger) is needed:
\(C(RTR^\dagger)C^\dagger=D\) implies
\((CR)T(CR)^\dagger=D\). An extra output permutation would only relabel
diagonal entries.

## Affine-router limitation and cost scope

For the two orientations of any unordered pair \(\{a,b\}\), their difference
in the original four label bits is \((d,d)\), where (d=a\oplus b\).
An invertible affine binary relabeling maps orientation differences by its
linear part. The three nonzero differences obey
\[
v_{11}=v_{10}+v_{01},\qquad v_d=(d,d).
\]
If every pair's two orientations were to differ in one bit after the
relabeling, the images of (v_{10}) and (v_{01}) would have to be distinct
unit vectors (they cannot coincide under an invertible map). Their sum has
Hamming weight two, so the image of (v_{11}) cannot be a unit vector.
Therefore no fixed affine CNOT/X relabeling makes even each pair separately
adjacent in one bit, much less assigns all six pairs a common orientation
bit. The proposed difference router is still a valid and cheap affine
relabeling; it leaves the (d=11) orientation as a two-bit flip.

This does not obstruct the controlled (C) above. It only shows that an
additional d-dependent orientation router, or a direct synthesis of the
sector-controlled (Q_d^\dagger), is needed if the chosen architecture
requires one common one-bit orientation coordinate. The abstract direct-sum
matrix is not yet a native gate factorization. The Bell-decoder prefix
contains two native F factors, whose generic compilation budget is four CX,
but at the exact decoder angles each F-plus-local block can be factored as
one CX; thus the exact prefix realization costs two CX. Adding the two-CX R
leaves eight CX in a 12-CX total budget when this exact prefix simplification
is used (or six CX under the generic compilation allocation). The entangler
cost of controlled C has not been derived here. No claim of a 12-CX
diagonalizer or obstruction to the fixed family follows.

## Common-phase route checked by hand

As a possible simplification of the controlled block basis, let
B = CS_A^dagger = diag(1, 1, 1, -i) on A. Applying B T_d B^dagger to |a>
gives these sector actions:

- d=00: CZ_A.
- d=10: i^(a1) X0.
- d=01: i^(a0) X1.
- d=11: i on even a parity and 1 on odd a parity, times X0 X1.

The d=10 and d=01 operators are diagonalized by H0 and H1, respectively.
The d=11 operator is diagonalized by the Bell decoder H0 CX(0->1): CX maps
X0 X1 to X0 and the parity phase to a diagonal function of the second bit,
after which H0 diagonalizes X0.

There is a conditional simplification. Apply P=CX(0->1) only when d1=1.
In d=01, both i^(a0) and X1 are unchanged by this conjugation. In d=11,
the parity phase becomes i^(1-a1) and X0 X1 becomes X0, so H0 now suffices.
Therefore after B and P, apply H0 when d0=1, and H1 when d0=0,d1=1.
Together with no mixing for d=00, these branch conditions give an exact
algebraic diagonalization when composed after R. This makes the desired
control structure explicit, but not its allowed-cost implementation: P is a
three-wire controlled-CX operation, and the conditional Hadamards and B
also need native/compiled decompositions. No cost within the remaining
budget is established.

## Handoffs

The exact block formulas and affine limitation are suitable inputs to a
separate constructive synthesis. The immediate synthesis target is the
explicit conditional route above (or equivalently C=direct-sum_d Q_d^dagger),
not a dense arbitrary eigenbasis and not the more general unordered-pair
code router. Any proposed cost must include control on d. The exact Bell
prefix (2 CX) plus R (2 CX) leaves 8 CX; counting the prefix at its generic
F compilation allocation (4 CX) plus R leaves 6 CX. Exact certification and
complete native/compiled gate-list checks remain necessary for any candidate.
