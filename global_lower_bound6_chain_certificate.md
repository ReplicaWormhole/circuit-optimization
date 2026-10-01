# A three-wire spectral-path obstruction for five-CNOT schedules

**Historical intermediate certificate.** The complete schedule proof in
`global_lower_bound6_complete_filter.md` establishes the current global
interval \(6\leq C_{\min}\leq14\). This certificate remains one part of it.

The global bound remains (5\leq C_{\min}\leq14). This note excludes one
further subclass of five-CNOT schedules after the two-terminal-axis filter
in `global_lower_bound6_two_terminal_axes.md`. It permits arbitrary local
gates, CNOT orientations and pairs, eigenvalue ordering, and eigenbases
within degenerate shift eigenspaces.

## Selected axis and its exact linear span

Suppose a degree-two input wire (j) has its second CNOT in position four
(counting chronologically from one), on pair (j,k). The fifth CNOT joins
(k,\ell), where (\ell\ne j). Choose the input Pauli axis on (j) that
commutes through its first incident CNOT. It remains local on (j) until
the fourth CNOT. Conjugation by that CNOT leaves its image in the span of
\(\sigma_\mu^{(j)}\) and
\(\sigma_\mu^{(j)}\sigma_\nu^{(k)}\), with
\(\mu,\nu\in\{X,Y,Z\}\): a CNOT never removes the nonidentity factor
from its original wire. The final CNOT sends a one-body Pauli on (k)
into the span of one-body terms on (k) and two-body terms on (k,\ell).
All two-body terms have the same Pauli axis on (\ell), even when the
CNOT direction is reversed. A collective one-qubit rotation, which
commutes with (V_4), makes that axis (Z_\ell).

For the physical chain ((j,k,\ell)=(0,1,2)), every selected final image
therefore lies in the **overapproximate** real space

\[
\mathcal L_{012}=\operatorname{span}_{\mathbb R}
\{\sigma_\mu^{(0)},\;\sigma_\mu^{(0)}\sigma_\nu^{(1)},\;
\sigma_\mu^{(0)}\sigma_\nu^{(1)}Z_2:
\mu,\nu\in\{X,Y,Z\}\},
\]

of dimension (3+9+9=21). This includes all local-gate choices and both
CNOT directions. It is larger than the actual two-CNOT-image family.

## Exact finite obstruction

For (M_a=4P_a=\sum_{r=0}^3i^{-ar}V_4^r), define

\[
F(A)=\sum_{a,b,c\;\mathrm{distinct}}
\lVert M_a A M_b A M_c\rVert_F^2.
\]

The spectral-path identity requires (F(A)=0) for every physical image
of an input one-qubit Pauli axis. Expand (A=\sum_{i=0}^{20}x_iB_i) in
the displayed Pauli basis, and let (m(x)) list the 231 monomials
(x_ix_j) for (i\leq j). The script
`global_lower_bound6_chain_gram.py` constructs the exact integer Gram
matrix (G) such that (F(A)=m(x)^\mathsf T Gm(x)). It is saved in
`global_lower_bound6_chain_gram_result.json` (ledger run 199). The rank
over (\mathbb Q) is 207, so matrix rank alone does not prove exclusion.

`global_lower_bound6_chain_certificate.py` makes the rank-one check with
exact rational arithmetic. If (F(A)=0), positivity of the Gram form
implies (Gm(x)=0). Exact nullspaces first force six (x_i^2) coordinates
to vanish, then two more. With the remaining 13 real coefficients, add
the normalization (\sum x_i^2=1) to the quadratic equations (Gm(x)=0).
The exact Gröbner basis over (\mathbb Q) is \(\{1\}\), so no nonzero real
(A\in\mathcal L_{012}) satisfies the spectral-path identity. The stored
result lists the eliminated coefficient indices and exact ranks at each
stage. The Gram entries originate in Gaussian-integer products. Each
quadratic image entry has magnitude at most 128 and each Gram entry is
bounded by (24\cdot256\cdot128^2<2^{27}), making the integer-valued
complex128 operations exact; the nullspace and Gröbner steps use SymPy
rationals.

Rotations and reflections of the physical four-cycle map every chain
with (j,k) adjacent and (j,\ell) opposite to (0,1,2). A reflection
sends (V_4) to (V_4^{-1}), merely permuting the spectral labels in
(F). Hence the obstruction applies to every such selected chain.

## Exhaustive schedule consequence

Among the 256 schedules remaining after run 195, exactly 48 have a
degree-two wire with this penultimate second incidence and final chain:
24 have five simple edges and 24 have one doubled edge. The script
enumerates all (6^5=7776) ordered unoriented pair schedules and leaves
208 schedules:

| Degree sequence and edge multiplicities | Remaining schedules |
|---|---:|
| (3322), five simple edges | 104 |
| (3322), one doubled edge and three simple edges | 64 |
| (4222), one doubled edge and three simple edges | 40 |

The remaining schedules are **not** validated diagonalizers. A next
calculation could apply the same exact rank-one test to a selected chain
where (j,k) and (j,\ell) are both physically adjacent, or to a
selected wire whose second CNOT occurs before position four. No global
six-CNOT lower bound is claimed.
