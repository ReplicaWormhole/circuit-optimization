# A same-pair output diagonal gauge cannot save a terminal CNOT

The exact 14-CNOT diagonalizer ends with two disjoint blocks
\(F_{01}\otimes F_{23}\), where

\[
F=\exp\bigl(i\pi(XX+YY)/8\bigr).
\]

An arbitrary diagonal two-qubit output gauge on either terminal pair is,
up to a global phase and local \(Z\) rotations,
\(G(c)=\exp(icZZ)\) for real \(c\). Every such gauge preserves
diagonalization, so one possible route to 13 CNOTs would replace a
two-CNOT \(F\) block by a one-CNOT circuit for \(G(c)F\).

This route is impossible. The Bell states simultaneously diagonalize
\(XX,YY,ZZ\). Since the three Pauli products commute, the eigenphases of
\(G(c)F\) in Bell order \((\Phi_+,\Phi_-,\Psi_+,\Psi_-)\) are

\[
c,\quad c,\quad \pi/4-c,\quad-\pi/4-c.
\]

In the magic basis, the local-equivalence invariant
\(m(U)=(Q^\dagger UQ)^\mathsf T(Q^\dagger UQ)\) consequently has
eigenvalues

\[
e^{2ic},\quad e^{2ic},\quad i e^{-2ic},\quad-i e^{-2ic}.
\]

The final two entries are distinct for every real \(c\). Hence this
multiset can never have the two-double-root pattern of CNOT, even after a
common global phase. It also cannot have the four-equal-root pattern of a
local gate. Local gates on either side preserve these multiplicities.
Therefore \(G(c)F\) needs at least two CNOTs for every \(c\).

This is an exact **fixed-pair** obstruction. A four-qubit diagonal output
gauge coupling the two terminal pairs, or a jointly resynthesized earlier
boundary, is outside its scope.
