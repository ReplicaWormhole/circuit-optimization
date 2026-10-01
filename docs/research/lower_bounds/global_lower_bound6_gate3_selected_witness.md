# Gate-three selected axes: exact surviving witnesses

**Historical intermediate witnesses.** The later half-shift trace
obstruction in `global_lower_bound6_complete_filter.md` rules out these
remaining cases and proves the global six-CNOT lower bound.

After the terminal and two penultimate-chain filters, 88 five-CNOT pair
schedules remain unexcluded. Forty-eight have an opposite-pair selected
axis at gate four, classified in
`global_lower_bound6_opposite_chain_certificate.md`. The other forty have
two degree-two wires whose second incidences occur at gates three and
five. An exhaustive scan in
`global_lower_bound6_gate3_selected_witness.py` separates these forty
into eight *adjacent-star* schedules and thirty-two *opposite-chain*
schedules. These are chronological unoriented-pair counts, with all
CNOT directions and one-qubit gates still free.

For the degree-two wire (j) finishing at gate three, choose an input
Pauli axis that commutes through its first CNOT. Its image remains
nonidentity on (j). In the adjacent-star representative
((j,k,\ell,m)=(2,3,0,1)), the selected axis next meets (j,k), then
(k,\ell), then (k,m). Arbitrary one-qubit gates and CNOT directions
put its image inside the overapproximate Pauli span with exact support
patterns

\[
\{j\},\quad\{j,k\},\quad\{j,k,\ell\},\quad
\{j,k,m\},\quad\{j,k,\ell,m\}.
\]

The dimensions of these mutually orthogonal Pauli-string sectors are
(3,9,27,27,81), totaling 147. In the opposite-chain representative
((j,k,\ell,m)=(2,0,1,3)), the later edges are (k,\ell) then
(\ell,m). The corresponding support patterns are
(\{j\},\{j,k\},\{j,k,\ell\},\{j,k,\ell,m\}\), of total dimension
(3+9+27+81=120). These spans are necessary overapproximations; they
do not impose the nonlinear restrictions of a particular CNOT gate.

Both spans contain the full-support Pauli product
(Z_0Z_1Z_2Z_3). The script goes further and supplies one **exact
five-CNOT circuit in each shape** for which

\[
U Z_2 U^\dagger=Z_0Z_1Z_2Z_3,
\qquad
U X_t U^\dagger=X_1X_3,
\]

with (t=1) in the adjacent-star representative and (t=3) in the
opposite-chain representative. The chronological CNOT lists, in
control-target order, are:

| Shape | CNOTs |
|---|---|
| Adjacent star | (0\to1,\ 2\to0,\ 3\to2,\ 0\to3,\ 1\to3) |
| Opposite chain | (0\to1,\ 2\to3,\ 0\to2,\ 1\to0,\ 3\to1) |

The images satisfy the exact spectral-path identity. The first commutes
with (V_4), as required for the diagonal input (Z_2). The second
commutes with (V_4^2), as required for an opposite-pair selected axis.
Because the diagonal input (Z_2) has full physical support, the
two-bit multiplicity argument does not apply. All equalities were
checked by integer matrix multiplication; the spectral projector
numerators use exactly represented Gaussian integers. The first
representative has 14 nonzero off-diagonal entries in (U^\dagger V_4U),
and so does the second. **Neither circuit diagonalizes the shift.**

This bounded test finds a rank-one spectral-path solution in each
gate-three support span and excludes **none** of the forty schedules.
It establishes the limit of this selected-axis test, not the viability
of a complete five-CNOT diagonalizer. A stronger invariant must use
additional input axes or the full commuting diagonal algebra.
