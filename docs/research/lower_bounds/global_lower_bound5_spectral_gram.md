# Global five-CNOT lower bound for diagonalizing the four-qubit shift

**Historical intermediate result.** The later complete five-CNOT schedule
proof in `global_lower_bound6_complete_filter.md` establishes the current
global interval \(6\leq C_{\min}\leq14\). The five-CNOT theorem below remains
a necessary step in that proof.

Let \(V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle\), and let
\(U^\dagger V_4 U=D\) be diagonal. An ancilla-free circuit for \(U\) using
arbitrary one-qubit gates and CNOTs between arbitrary pairs needs at least
**five CNOTs**. The exact 14-CNOT construction gives the global interval
\(5\leq C_{\min}\leq14\). No output eigenvalue ordering or eigenbasis inside
a degenerate sector is fixed.

The spectral-path identity in
`global_lower_bound4_spectral_path.md` says that, for every input one-qubit
Pauli axis \(p\), its physical image \(A=UpU^\dagger\) satisfies

\[
P_a A P_b A P_c=0\quad\text{for all pairwise distinct }a,b,c\in\{0,1,2,3\},
\]

where (P_k) is the (i^k) spectral projector of (V_4). This allows arbitrary
one-qubit gates and arbitrary eigenbases inside degenerate sectors. It is the
four-qubit analogue of the spectral-path constraint in Appendix C.1, Lemma 11
of the supplied *Exact ancilla-free diagonalization of qubit cycles* PDF.

## Exact finite certificate

Any physical one-CNOT image of a one-qubit Pauli axis, supported on a pair
(a,b), has the form

\[
A=\alpha\,\sigma(\mathbf u)_a+\beta\,\sigma(\mathbf v)_a
                       \sigma(\mathbf w)_b,
\]

up to swapping (a,b). This follows by conjugating a local Pauli axis through
a maximal Ising interaction locally equivalent to CNOT. A collective rotation
(R^{\otimes4}), which commutes with (V_4), puts
(\sigma(\mathbf w)_b=Z_b). Thus (A) belongs to the real span

\[
\mathcal L_{ab}=\operatorname{span}_{\mathbb R}
  \{X_a,Y_a,Z_a,X_aZ_b,Y_aZ_b,Z_aZ_b\}.
\]

Put (M_k=4P_k=\sum_{r=0}^{3}i^{-kr}V_4^r), and define the nonnegative
quartic form

\[
F(A)=\sum_{a,b,c\,\mathrm{distinct}}\|M_a A M_b A M_c\|_F^2.
\]

For coefficients (x_0,\ldots,x_5) in the displayed basis, let (m(x))
be the vector of the 21 monomials (x_i x_j) with (i\leq j). Exact
Gaussian-integer matrix multiplication yields (F(A)=m(x)^\mathsf T Gm(x)).
For the adjacent pair (01), `spectral_path4_adjacent_gram.py` saves the
integer (21\times21) matrix (G) in
`spectral_path4_adjacent_gram_result.json`. Its exact rank is 21. Since it is
a Gram matrix, it is positive definite; hence (F(A)>0) for every nonzero
(A\in\mathcal L_{01}). The other adjacent pairs and endpoint choices follow
by cyclic translation and a reflection of the four-wire cycle: reflection
conjugates (V_4) to (V_4^{-1}), merely permuting the spectral labels
summed in (F). The exact eigenvalues provide a second compact
certificate: (12288) (twice), (16384), (20480) (three times), (24576),
(32768) (seven times), (40960) (five times), and
(81920\pm32768\sqrt2). This is ledger run 185. Although the script uses
complex128 arrays, these operations are exact Gaussian-integer arithmetic:
each row of every \(M_k\) has absolute row sum at most four, each Pauli
matrix is monomial, and each quadratic image has absolute row sum at most
128. Every Gram entry is therefore bounded by
\(24\cdot16^2\cdot128^2<2^{27}\), far below the integer-exact range of
IEEE binary64.

For the opposite pair (02), the analogous Gram matrix has rank 20 and
one-dimensional nullspace, spanned by the monomial (x_5^2), where
(x_5) is the coefficient of (Z_0Z_2). Consequently a nonzero one-CNOT
image on an opposite pair obeying all spectral-path identities must be
proportional to \(\sigma(\mathbf n)_0\sigma(\mathbf n)_2\) for a single
Bloch axis \(\mathbf n\). This follows by undoing the collective rotation.
The matrix and nullspace are in `spectral_path4_opposite_gram_result.json`
(ledger run 187). Same-axis antipodal products do survive because they
commute with (V_4^2).

## Excluding four CNOTs

The incidence argument in `global_lower_bound4_spectral_path.md` forces
every wire of a hypothetical four-CNOT diagonalizer to meet exactly two
CNOTs. The supplied paper's Section 7, Lemma 10 and proof of Theorem 2
exclude a disconnected interaction graph, leaving a four-cycle.

Consider either endpoint of the chronologically last CNOT of (U). Choose
an input Pauli axis that commutes through its one earlier incident CNOT.
Its physical image encounters only the last CNOT after that and is therefore
a one-CNOT image supported on the last pair. The adjacent-pair certificate
forbids the last pair from being adjacent under the physical shift. Hence
**the final CNOT of a four-CNOT candidate must join opposite physical wires**:
\(02\) or \(13\). Its two selected endpoint observables must be same-axis
antipodal products. Each such product commutes with \(V_4^2\), since the
half-shift exchanges the opposite wires.

Write \(p_j\) for the selected input axis on each endpoint \(j\) of the
last CNOT, and \(A_j=Up_jU^\dagger\). The physical support of \(A_j\) is
contained in the two-wire last pair. Therefore \(p_j\ne\pm Z_j\):
otherwise \(A_j\) would commute with \(V_4\), because \(D\) is diagonal,
and Lemma 10 would force this non-scalar \(A_j\) to have full four-wire
support. Each \(p_j\) consequently has a nonzero off-diagonal matrix
element between the two computational states of its input wire.

Since \(A_j\) commutes with \(V_4^2\), its preimage \(p_j\) commutes with
\(D^2\). The latter is diagonal. The off-diagonal entry of \(p_j\) makes
the two diagonal entries of \(D^2\) equal across *every* edge that flips
bit \(j\). This holds for the two distinct endpoint bits of the last
CNOT. Thus every eigenvalue multiplicity of \(D^2\) is divisible by four.
But \(V_4^2\) swaps wires \(0\leftrightarrow2\) and
\(1\leftrightarrow3\), so \(\operatorname{tr}V_4^2=4\). Its \(+1\) and
\(-1\) eigenspaces have dimensions \((16+4)/2=10\) and \((16-4)/2=6\),
neither divisible by four. This contradiction rules out four CNOTs.
The preceding incidence argument already rules out fewer than four.

## Verification scope

- `python3 archive/legacy_root/spectral_path4_adjacent_gram.py` (run 185): exact integer Gram rank 21.
- `python3 archive/legacy_root/spectral_path4_opposite_gram.py` (run 187): exact integer Gram rank 20,
  nullspace spanned by (x_5^2).
- `python3 archive/legacy_root/spectral_path4_pair_pauli_scan.py` (run 182): Pauli-product diagnostic
  only; all adjacent-pair Pauli products fail, while (XX,YY,ZZ) survive on
  opposite pairs.

The exact-rank arguments are analytic once the finite matrices are specified.
The proof applies to arbitrary local gates, connectivity, output ordering,
and eigenbasis choices. It does not show that the 14-CNOT construction is
optimal.
