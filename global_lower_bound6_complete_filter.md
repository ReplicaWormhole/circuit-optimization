# Global six-CNOT lower bound for diagonalizing the four-qubit shift

Let $V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle$. Let $U$ be an
ancilla-free circuit with arbitrary one-qubit unitaries and CNOTs between
arbitrary physical pairs, and suppose $U^\dagger V_4U=D$ is diagonal.
Every such circuit contains **at least six CNOTs**. The exact 14-CNOT
construction gives $6\leq C_{\min}\leq14$. The claim fixes neither
the order of the eigenvalues of $D$ nor bases within degenerate
eigenspaces.

## A half-shift trace obstruction

Let $p_j=\mathbf n\cdot\boldsymbol\sigma_j$ be an input one-qubit
Pauli axis and $A=Up_jU^\dagger$. Suppose $A$ has proper physical
support and commutes with $V_4^2$. Then $p_j$ is not diagonal: if it
were $\pm Z_j$, $A$ would commute with $V_4$ because $D$ is
diagonal, contradicting the full-support lemma. Indeed, the minimal
physical support of an operator commuting with $V_4$ is invariant
under cyclic translation of the four wires; the only invariant subsets
are empty and all four. A non-scalar operator must therefore have full
support. Hence $p_j$ has a nonzero off-diagonal entry. From
$[p_j,D^2]=0$, the diagonal entries of $D^2$ agree for every pair
of input basis strings differing only in bit $j$. Pairing these strings
in the trace gives

\[
\operatorname{tr}(V_4^2A)=\operatorname{tr}(D^2p_j)=0.\tag{1}
\]

But if $a,b$ are opposite physical wires, $V_4^2$ swaps $a,b$ and
also swaps the complementary opposite pair. For every unit Bloch axis
$\mathbf n$, the swap trace identity yields

\[
\operatorname{tr}\!\left(V_4^2
\sigma(\mathbf n)_a\sigma(\mathbf n)_b\right)
=\operatorname{tr}(\mathrm{SWAP}_{ab}
\sigma(\mathbf n)\otimes\sigma(\mathbf n))
\operatorname{tr}(\mathrm{SWAP}_{cd})=2\cdot2=4.\tag{2}
\]

Thus **no input one-qubit Pauli axis can have a proper-support image
equal to $\pm\sigma(\mathbf n)_a\sigma(\mathbf n)_b$ on opposite
physical wires**. `global_lower_bound6_complete_filter.py` checks the
underlying $3\times3$ exact trace tensor for both opposite pairs
using 16-by-16 SymPy matrices: it is $4I_3$ in the $X,Y,Z$ basis.

## Finite schedule closure

The incidence argument in `global_lower_bound4_spectral_path.md`
proves that every input wire must meet at least two CNOTs. The
full-support lemma requires the forward light cone of each input wire
to cover all four physical wires, because $UZ_jU^\dagger$ commutes
with $V_4$. We now enumerate all $6^5=7776$ chronological sequences
of **unoriented** CNOT pairs; every CNOT direction remains allowed.

Call a degree-two wire's second incident CNOT *terminal* when no later
CNOT touches either endpoint of that pair. A commuting input axis passes
through its first incident CNOT, then has a one-CNOT image on the
terminal pair. The exact spectral-path Gram classification in
`global_lower_bound5_spectral_gram.md` forbids an adjacent terminal
pair. On an opposite terminal pair, that classification forces the
image to be a same-axis product, now forbidden by (1)-(2). **A single
degree-two endpoint of a terminal edge suffices for exclusion.**

Every schedule surviving those terminal tests has a degree-two wire
$j$ whose second incident CNOT is fourth, on pair $j,k$, and whose
fifth CNOT touches $k$ but not $j$, on pair $k,\ell$. There are
three physical chain shapes, each represented 48 times:

1. $j,k$ adjacent and $j,\ell$ opposite. The exact 21-dimensional
   overapproximate Pauli-image Gram and rank-one Gröbner certificate in
   `global_lower_bound6_chain_certificate.md` show that **no nonzero**
   image in this class satisfies the spectral-path identity.
2. Both $j,k$ and $j,\ell$ adjacent. The separate exact
   21-dimensional rank-one certificate in
   `global_lower_bound6_adjacent_v_certificate.md` again has **no
   nonzero** spectral-path solution.
3. $j,k$ opposite. The exact rank-one classification in
   `global_lower_bound6_opposite_chain_certificate.md` forces every
   real spectral-path image in the 21-dimensional overapproximation to
   have support only on opposite wires $j,k$, with a symmetric real
   $3\times3$ Pauli coefficient matrix $M$. Because the actual image
   is a Hermitian involution, collectively rotate to diagonalize $M$:
   $A=\lambda_X X_jX_k+\lambda_Y Y_jY_k+\lambda_Z Z_jZ_k$.
   Squaring gives

   \[
   A^2=(\lambda_X^2+\lambda_Y^2+\lambda_Z^2)I
   -2\lambda_Y\lambda_Z X_jX_k
   -2\lambda_X\lambda_Z Y_jY_k
   -2\lambda_X\lambda_Y Z_jZ_k.
   \]

   Since $A^2=I$, at most one coefficient is nonzero, and it is
   $\pm1$. Thus $A=\pm\sigma(\mathbf n)_j\sigma(\mathbf n)_k$,
   contradicted by (1)-(2).

The exclusive exhaustive counts are:

| First necessary obstruction met | Schedules |
|---|---:|
| Some wire has fewer than two CNOT incidences | 5,196 |
| A commuting input $Z_j$ cannot reach all four wires | 1,572 |
| Terminal adjacent edge with a degree-two endpoint | 576 |
| Terminal opposite edge with a degree-two endpoint | 288 |
| Penultimate adjacent chain with opposite endpoints | 48 |
| Penultimate adjacent V-shaped chain | 48 |
| Penultimate opposite-first chain | 48 |
| **Total** | **7,776** |

`global_lower_bound6_complete_filter.py` reproduces this partition and
asserts that no schedule remains (ledger run 223). The 16-by-16 trace
calculation uses exact SymPy arithmetic. The two adjacent-chain Gram
matrices have exact Gaussian-integer entries; their rank-one Gröbner
calculations use rational arithmetic and include the normalization
$\sum_i x_i^2=1$, automatic for a normalized Pauli image expanded in
the Pauli basis. The prior five-CNOT theorem already excludes four or
fewer CNOTs, so this exhaustive five-CNOT exclusion proves the stated
global lower bound.

## Scope and historical records

This is a lower bound for diagonalizing $V_4$, not for the stronger
Schur transform. The 14-CNOT upper construction was independently
checked earlier; it remains an upper bound rather than an optimality
proof. The experiment ledger and intermediate notes remain useful
records of how the finite certificates were found. Statements there
that a *single* terminal opposite axis survives, or that the current
global lower bound remains five, are superseded by the trace lemma and
the complete schedule closure above.
