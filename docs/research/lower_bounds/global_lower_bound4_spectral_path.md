# Global four-CNOT lower bound for diagonalizing the four-qubit shift

Let (V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle). We allow arbitrary
one-qubit gates, all CNOT pairs, every output ordering, and every orthonormal
basis and phase choice within degenerate shift eigenspaces. If an ancilla-free
circuit (T) obeys (T V_4 T^\dagger=D) with (D) diagonal, then it contains
**at least four CNOTs**. This is an intermediate bound; the stronger
five-CNOT theorem in `global_lower_bound5_spectral_gram.md` uses the same
spectral-path identity. The subsequent complete five-CNOT exclusion in
`global_lower_bound6_complete_filter.md` gives the current interval
(6\leq C_{\min}\leq14); optimality remains open.

## Spectral-path obstruction

Use the inverse circuit (U=T^\dagger), so (U^\dagger V_4U=D). Let (P_k)
be the spectral projector of (V_4) for eigenvalue (i^k). For every
one-qubit Pauli axis (p=\vec v\cdot\vec\sigma) acting on an *input* wire of
(U), put (A=UpU^\dagger). Then

\[
  P_a A P_b A P_c=0\qquad(a,b,c\text{ pairwise distinct}).
\]

Indeed, in the input computational basis each spectral projector of (D)
is diagonal and (p) only connects a bit string to itself or to the string
with its one selected bit flipped. To traverse two distinct eigenvalue
labels, both applications of (p) must flip the bit; after the second flip
the string has returned to its original label. Thus a path through three
distinct labels vanishes. Conjugation by (U) gives the displayed identity.

No nonzero physical one-qubit Pauli axis satisfies this identity. To see it
for (Z_0), restrict to the one-excitation subspace. Write (e_j) for the
string with its sole 1 at wire (j), and choose Fourier eigenvectors

\[
  \phi_k=\frac12\sum_{j=0}^3 i^{-kj}e_j,\qquad k=0,1,2,3.
\]

The shift moves (e_j\mapsto e_{j+1\bmod4}), so (V_4\phi_k=i^k\phi_k).
On this invariant subspace (Z_0=I-2|e_0\rangle\langle e_0|), whence
(\langle\phi_a|Z_0|\phi_b\rangle=-1/2) for (a\ne b). Since each shift
eigenspace is one-dimensional *within this weight sector* and (Z_0)
preserves weight,

\[
  \langle\phi_0|P_0 Z_0P_1 Z_0P_2|\phi_2\rangle=\tfrac14\ne0.
\]

Any one-qubit Pauli axis on any physical wire is conjugate to (Z_0) by a
collective one-qubit rotation and a power of (V_4). Both commute with the
spectral projectors, so the same nonvanishing follows for every axis.

## Counting CNOT incidences

Suppose an input wire of (U) meets no CNOT. Any Pauli axis on that wire
remains physically one-body under conjugation by (U), contradicting the
spectral-path obstruction. If it meets exactly one CNOT, choose its input
Pauli axis so that immediately before that CNOT it is (Z) on a control or
(X) on a target. This is possible because all earlier gates on that wire
are one-qubit gates. The chosen axis commutes through the CNOT and all later
CNOTs are on other wires; hence its final image is again one-body, the same
contradiction. Thus every one of the four wires has at least two CNOT
incidences. Each CNOT contributes two incidences, proving (C\geq4).

With exactly four CNOTs, the interaction multigraph has degree two at each
vertex. Lemma 10 and the proof of Theorem 2 in Section 7 of the supplied
*Exact ancilla-free diagonalization of qubit cycles* paper give the
full-support commuting-observable argument that excludes a disconnected
multigraph, leaving a four-cycle. The accompanying exact
finite script checks a weaker necessary condition: eight of the 24 edge
schedules of a four-cycle allow every output wire's backward light cone to
reach all four wires. Therefore this condition alone does not exclude four
CNOTs. The stronger spectral Gram argument in
`global_lower_bound5_spectral_gram.md` excludes those remaining schedules;
raising the bound above five remains open.

Run `python3 archive/legacy_root/global_lower_bound4_spectral_path.py` for the rational Fourier
matrix-element check and the four-cycle schedule count. The mathematical
proof above does not rely on floating-point computation.
