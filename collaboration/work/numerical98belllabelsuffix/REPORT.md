# Bell-label exchange eigenbasis and reordered-suffix hypothesis

Board 109. Structural derivation only. I read the Bell-prefix report and the
saved run402 result; I performed no numerical or symbolic matrix evaluation,
search, optimization, or fitting.

## Exact exchange blocks

Write the two Bell-register labels as two-bit strings (a=(a_0,a_1)) on
(A=(q_0,q_2)) and (b=(b_0,b_1)) on (B=(q_1,q_3)). From

\[
T=(CZ_A\otimes I_B)\,SWAP_{AB},
\]

the action on the Bell-label basis is

\[
T|a,b\rangle=\sigma_b|b,a\rangle,
\qquad \sigma_x=(-1)^{x_0x_1}.
\]

Thus (|a,a\rangle) is already an eigenvector with eigenvalue \(\sigma_a\).
For a distinct unordered pair \(\{a,b\}\), the invariant block in ordered
basis \((|a,b\rangle,|b,a\rangle)\) is

\[
T_{ab}=\begin{pmatrix}0&\sigma_a\\\sigma_b&0\end{pmatrix},
\qquad T_{ab}^2=\sigma_a\sigma_b I.
\]

When \(\sigma_a=\sigma_b=s\), its eigenvectors are

\[
|\psi_{ab}^{\pm}\rangle
=\frac{|a,b\rangle\pm|b,a\rangle}{\sqrt2},
\qquad \lambda_\pm=\pm s.
\]

When the signs differ, for each \(\lambda\in\{+i,-i\}\),

\[
|\psi_{ab}^{\lambda}\rangle
=\frac{|a,b\rangle+(\lambda/\sigma_a)|b,a\rangle}{\sqrt2}.
\]

There are three positive labels \(00,01,10\) and one negative label \(11\).
The resulting spectrum is \(+1\) with multiplicity 6 (three positive
diagonal states and three symmetric positive-positive blocks), \(-1\) with
multiplicity 4 (the \(|11,11\rangle\) diagonal state and three
antisymmetric positive-positive blocks), and \(+i,-i\) each with multiplicity
3 (the three positive/negative blocks).

Let (Q) be the unitary whose columns are these eigenvectors, arranged inside
each singleton or two-dimensional invariant block. Then

\[
W=Q^\dagger,\qquad WTW^\dagger=D.
\]

On each positive-positive two-state block, (Q=H), so (W=H). On a block
ordered with positive-sign label first and (11) second, (Q=SH), where
\(S=\operatorname{diag}(1,i)\), so (W=HS^\dagger). Singleton blocks need no
mixing. These are six independent weighted exchange rotations: three (H)
blocks and three phase-weighted (HS^\dagger) blocks. They recover the
\((6,4,3,3)\) root multiplicities without making an arbitrary dense-basis
assumption.

## Constructive routing target and architecture hypothesis

A constructive way to realize this block unitary in a general four-qubit
circuit is to define a reversible label router (R). For each distinct
unordered pair, map its two orientations into one dedicated adjacent
two-state code block; map the four diagonal labels into four remaining
singleton codes. For example, assign the six unordered pairs to the twelve
codes 0 through 11 as three-bit block index plus one orientation bit, and
assign the four diagonal states to codes 12 through 15. In that routed basis,
apply (H) to the orientation bit for the three equal-sign blocks and
\(HS^\dagger\) for the three opposite-sign blocks, controlled on the block
index, then undo (R). This gives the explicit synthesis target

\[
W=R^\dagger\,C_{\rm weighted\ exchange}\,R.
\]

This is an exact algorithmic construction of (W), but no low-cost
decomposition of (R) or the controlled block operation into the reserved
suffix has been established. The reserved suffix after the Bell-decoder prefix
is

`CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.

In the (A/B) label bits, the first three CXs form a parity fan-in at (q_1):
they update \(b_0\mapsto b_0\oplus a_0\oplus a_1\oplus b_1\). `CX12` then
feeds that parity into (a_1). The final native gates `F01` and `F23` act on
the cross-register bit pairs \((a_0,b_0)\) and \((a_1,b_1)\). This makes a
weighted-exchange routing hypothesis concrete: use the fan-in and the free
interior rotations (A1,B1) to encode/control orientation and sign class,
then use the cross-register F gates and remaining local layers to implement
the block rotations and cleanup. The three equal-sign blocks require
unweighted exchange mixing; the three opposite-sign blocks require an
additional conditional \(S^\dagger\) phase. This is a guide for a new
constructive initializer or synthesis attempt, not a demonstrated assignment
of the six native suffix entanglers (four CX and two F gates), which cost
eight compiled CXs. The Bell prefix has two additional F gates, adding four
compiled CXs for a total compiled cost of twelve.

The run402 record reports a single bounded fit ending at loss 0.375 with an
invalid 12-CX compiled list. That plateau motivates changing the suffix
initialization from an identity-like residual start to the routed, weighted
exchange target. It supplies no proof that the reordered family is incapable
of diagonalization. A realistic next step is to synthesize the explicit router
and controlled block operation, calculate its actual native cost, and only
then decide whether the fixed suffix can host it. If it costs more than the
available four CX plus two F gates after the Bell prefix, retain it as a
structural comparator rather than claiming this family has an exact
12-CX realization. Any eventual numerical candidate still needs independent
full-list validation and exact certification.
