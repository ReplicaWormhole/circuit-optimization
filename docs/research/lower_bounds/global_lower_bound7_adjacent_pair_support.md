# General adjacent-pair obstruction and six-CNOT schedule reduction

This note extends the one-CNOT adjacent-pair spectral obstruction to
**every traceless Hermitian operator supported on an adjacent physical
pair**. It applies to arbitrary local one-qubit gates and both CNOT
directions. It does not establish a global seven-CNOT lower bound.

## Exact spectral-path certificate

Let \(P_a\) be the \(i^a\)-eigenspace projector of the four-qubit
cyclic shift \(V_4\), and \(M_a=4P_a\). For an image \(A=Up_jU^\dagger\)
of any input one-qubit Pauli axis in a diagonalizer, the spectral-path
identity requires

\[
F(A)=\sum_{a,b,c\,\mathrm{distinct}}
\|M_a A M_b A M_c\|_F^2=0.
\]

For the physical adjacent pair \((0,1)\), expand \(A\) in all 15
nonidentity two-wire Pauli strings, with real coefficients \(x_i\).
`global_lower_bound7_adjacent_pair_general_gram.py` constructs the
exact integer Gram matrix \(G\) on the 120 quadratic monomials
\(m_{ij}=x_ix_j\), such that \(F(A)=m^\mathsf TGm\). The 120-by-120
matrix has exact rational rank 111. Every Gram entry is an integer of
absolute value at most 81,920, below the conservative exact-binary64
bound \(2^{27}\) for this Gaussian-integer matrix computation.

Since \(G\) is a Gram matrix, \(F(A)=0\) implies \(Gm=0\). Its exact
nullspace forces the squares of all six one-body Pauli coefficients
to vanish. For a real Hermitian \(A\), those coefficients are zero.
`global_lower_bound7_adjacent_pair_general_certificate.py` substitutes
this into \(Gm=0\), adds \(\sum_i x_i^2=1\) for the remaining nine
two-body coefficients, and computes the Gröbner basis \(\{1\}\) over
\(\mathbb Q\). Thus **no nonzero traceless Hermitian adjacent-pair
operator satisfies the spectral-path identity**. Any nonzero real
solution could be rescaled to the normalization, so the certificate
is not restricted to Pauli involutions. Dihedral symmetry carries
the result to every adjacent physical pair; reflection permutes the
spectral labels summed in \(F\).

The Gram and corrected rank-one certificate are ledger runs 237 and
240. Run 239 was marked failed: its basis was correctly \(\{1\}\),
but a coefficient-domain comparison produced an incorrect summary
boolean. Run 240 uses `gb.contains(1)` and its result reports no
normalized nonzero solution.

## Orientation-independent support screen

For each input wire \(j\), select a Pauli axis commuting through its
first incident CNOT. Starting with potential support \(\{j\}\)
immediately after that CNOT, scan all later CNOTs. Whenever an edge
intersects the potential support, add both endpoints. This cone
contains the actual final support for every choice of one-qubit gates
and CNOT directions. If the final cone is one adjacent physical pair,
the selected nonzero image contradicts the theorem above.

`global_lower_bound7_degree3_tail_screen.py` first classifies the
3,576 schedules surviving the earlier six-CNOT filter. Of the 1,416
degree-3333 schedules, 1,056 have a degree-three wire whose selected
axis meets exactly two active CNOTs after its first incidence. The
existing two-CNOT-chain certificates do not automatically cover the
degree-three \(j\)-centered or repeated-pair shapes, so that inventory
alone excludes nothing (ledger run 236).

`global_lower_bound7_adjacent_pair_support_filter.py` then applies
the new general adjacent-pair theorem to **all** 3,576 earlier
survivors (ledger run 242):

| Degree sequence | Earlier survivors | Newly excluded | Remain |
|---|---:|---:|---:|
| 3333 | 1,416 | 144 | 1,272 |
| 4332 | 1,896 | 144 | 1,752 |
| 4422 | 144 | 0 | 144 |
| 5322 | 120 | 80 | 40 |
| **Total** | **3,576** | **368** | **3,208** |

These are ordered schedules of unoriented CNOT pairs. A remaining
schedule is **not** a verified diagonalizer. No claim is made about
the still-open 3,208 schedules; the established global interval
remains \(6\le C_{\min}\le14\).

For comparison, `global_lower_bound7_opposite_pair_general_gram.py`
computes the analogous exact Gram matrix on the full opposite-pair
Pauli space (ledger run 243). It has rank 75 among 120 quadratic
monomials, with no individual coefficient square forced to zero by
its nullspace. That calculation supplies no opposite-pair exclusion;
the 3,208 count does not use it.
