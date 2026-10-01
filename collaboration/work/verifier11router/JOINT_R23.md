# Joint router target and restricted synthesis obstruction

Board130; independent by-hand verification of algebraic_sol/root identities.
No numerical matrices, symbolic search or optimizer were evaluated.
Let x=q0,y=q2,u=q1,v=q3 and R=CX(y->v). The mergedBP M preserves x and
has M0=I,M1=CX(v->y) S_y CZ_yv. This reproduces x1v0 S and x1v1 XSdag.
For D=S_yCZ, its diagonal phase after R conjugation is
 i^y (-1)^(y*(v xor y)) = (-i)^y (-1)^(yv),
so R D R=Sdag_y CZ. Also R CX(v->y) R=SWAP_yv. Therefore
 R M1 R=SWAP_yv CZ_yv Sdag_y,
with exact relative phases; order is chronological Sdag,CZ,SWAP.
M0 remains I. This is a new controlled phase-weighted fermionic-swap target,
not literal deletion, and can mix v sectors. Its synthesis cost is unproved.

The final selector F has u1 I and u0 (Y_y for v0,H_y for v1). Conjugating
by R gives F'=R F R. Its u0 action, in |y,v> order, is
00->i11,11->-i00,01->(01+10)/sqrt2,10->(01-10)/sqrt2.
Thus its active pair operator T' has even block Y and odd block H. Both
blocks have determinant-1, so it is a parity-preserving matchgate target.
F' is I for u1. CH commutes with R. A potential joint11CX circuit needs
honest costs for both the transformed merged block and final selector;
a four-CX merged block alone does not supply the missing saving.

## Restricted three-CX conjugation route is impossible

Root proposed testing a pair diagonalizer K with one CX plus arbitrary
one-qubit gates, a middle spectator-controlled local Pauli (one CX), and
Kdag (one CX). This requires T'=K(P tensor I)Kdag or the reversed-wire
version for a one-qubit Pauli-axis involution P.

The exact Pauli expansion is
 T'=(X tensor Y+Y tensor X)/2
    +(X tensor X+Y tensor Y+Z tensor I-I tensor Z)/(2sqrt2).
The Y-even sign is confirmed by 00->i11 and11->-i00. In the local
operator bases(I,X,Y,Z), the X/Y coefficient submatrix has diagonal
1/(2sqrt2) and offdiagonal1/2, with determinant1/8-1/4=-1/8.
The independent Z tensor I and I tensor Z entries add two further ranks.
Therefore the coefficient matrix, equivalently the operator Schmidt rank
across y|v, has rank4.

Conjugating any single-qubit Pauli-axis operator by one CNOT has rank at
most2. For a target-local axis, X maps to I tensor X while Y,Z map to
Z tensor Y,Z; for a control-local axis, Z maps to Z tensor I while X,Y
map to X,Y tensor X. Arbitrary local gates before and after the CNOT
rotate these local operator factors and preserve their Schmidt rank.
Hence K(P tensor I)Kdag has rank<=2, contradicting rank4 of T'.

This excludes only the stated one-CX pair-diagonalizer/middle-local-Pauli
three-CX construction. It does not exclude every three-CX implementation
of the controlled matchgate, joint router/selector synthesis, or an11CX
V4 diagonalizer. No generic lower bound follows. Adversarial review checked
signs, target/control conjugation cases, rank versus coefficient determinant,
and this restricted scope. No unresolved issue found in the derivation.
