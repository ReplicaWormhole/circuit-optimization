# Cross-pair parity phases on the fixed terminal suffix

Let $S=F_{01}\otimes F_{23}$, where
$F=\exp(i\pi(XX+YY)/8)$. This is the final pair of blocks in the exact
14-CNOT circuit, and the existing synthesis uses four CNOTs. For each
nonempty parity mask $m$ meeting both pairs — nine masks in total —
consider the output diagonal gauge

\[
G_{m,k}=\exp(i k\pi Z_m/16),\qquad k=0,\ldots,15.
\]

The remaining $k$ values modulo 32 change only a global sign. Multiplying
the complete diagonalizer on the left by $G_{m,k}$ preserves its cycle
diagonalization. This note concerns the CNOT cost of the **fixed suffix**
$G_{m,k}S$, with arbitrary one-qubit gates available on both sides.

`algebraic13_crosspair_parity_gauge_rank.py` computes operator-Schmidt ranks
for the four one-wire cuts and three two-versus-two cuts. It specializes the
entries in both \(\mathbb F_{193}\) and \(\mathbb F_{257}\), each containing a
primitive 32nd root of unity. Nonzero minors after specialization certify
nonzero minors over the cyclotomic field. A separate complex-matrix rank
calculation matched all 144 modular rank profiles.

| Angles | Cases | Ranks for cuts `0`, `1`, `2`, `3`, `01`, `02`, `03` | Cut-graph lower bound |
|---|---:|---|---:|
| $k=0,8$ | 18 | `(4,4,4,4,1,16,16)` | 4 CNOTs |
| All other sampled $k$ | 126 | `(4,4,4,4,2,16,16)` | 5 CNOTs |

For a circuit with at most four CNOTs, rank 16 across each of `02|13` and
`03|12` requires every CNOT to cross both cuts. Only edges `01` and `23`
cross both. A circuit using only those edges is local across `01|23` and
therefore has rank 1 there. This contradicts the rank 2 certified for every
nontrivial sampled gauge. Thus these fixed gauged suffixes require at least
five CNOTs. At $k=0,8$, the gauge is a product of one-qubit Pauli gates
up to phase, so the original four-CNOT implementation remains sufficient.

This modular screen gives an exact rank obstruction for the sampled angles.
The symbolic calculation below extends it to every real angle for a single
mask. Neither calculation covers products of multiple parity gauges or
moving the suffix boundary into an earlier block. The five-CNOT lower
bound is not a five-CNOT synthesis or a lower bound on the full four-cycle
diagonalizer.

## All real angles for single masks

`algebraic13_crosspair_parity_symbolic.py` extends the five masks of weight
three or four to every real angle $c$. After removing the nonzero global
even-parity phase from $\exp(icZ_m)$, the odd-parity output rows carry
$t=e^{-2ic}$. For each of the cuts `02|13` and `03|12`, the exact
realignment determinant is

\[
\det R(GS)=
\begin{cases}
\dfrac{t^4(t^2-2)^2(2t^2-1)^2}{65536},
&m\in\{7,11,13,14\},\\[4pt]
-\dfrac{(2t^2-1)^4(8t^2-9)}{65536},&m=15.
\end{cases}
\]

All zeros have modulus different from one, while a physical real $c$
always gives $|t|=1$. Thus both balanced cuts have rank 16 for every
real $c$. Across `01|23`, $S$ is local and the gauge has Schmidt rank
two unless $c\in(\pi/2)\mathbb Z$; in those exceptional cases the gauge
is a product of one-qubit Pauli gates up to phase. The same four-edge
contradiction therefore gives a **five-CNOT lower bound on this fixed
suffix for every nontrivial real-angle weight-three or weight-four
single-mask gauge**. Direct complex determinants at 40 mask/cut/angle
points matched the symbolic factors to at most $3.13\times10^{-17}$.

The four two-body cross-pair masks $m\in\{5,6,9,10\}$ also have a simple
all-angle determinant. For both balanced cuts,

\[
\det R(GS)=\frac{t^8}{65536}.
\]

This is nonzero for every physical $|t|=1$. The same rank-two argument
across the terminal-pair cut therefore extends the five-CNOT fixed-suffix
lower bound to **all nine single cross-pair parity masks at every
nontrivial real angle**. At $c\in(\pi/2)\mathbb Z$, the gauge is local
and the cost is exactly four CNOTs. The two-body factors were independently
checked against 32 direct complex determinant evaluations, with maximum
absolute difference $4.53\times10^{-20}$. The exact factorization is
reproduced by algebraic13_crosspair_two_body_symbolic.py.

This all-angle conclusion does not cover products of multiple masks or
shifted suffix boundaries.

## Two coordinated masks on a finite grid

The successor screen in algebraic13_crosspair_two_mask_grid.py considers
all 36 unordered pairs of distinct cross-pair masks. Each angle numerator
is independently chosen from $\{-4,-2,-1,1,2,4\}$ in units of $\pi/16$,
giving 1,296 fixed-suffix gauges. Exact ranks over both finite fields
matched an independent complex rank calculation for every case.

All 1,296 cases have ranks 4 on each one-wire cut and 16 on both balanced
cross cuts. On the terminal-pair cut, 648 have rank 2 and 648 have rank 4.
The same four-edge contradiction therefore certifies a **five-CNOT lower
bound for each sampled two-mask fixed suffix**. None is three-CNOT
compatible, and there was no candidate to synthesize.

The finite grid does not prove the corresponding statement for arbitrary
pairs of real angles. A concrete next test is to factor the two balanced-cut
realignment determinants as bivariate polynomials in the two relative
phases and check for zeros on the two-dimensional unit torus. A changed
boundary reaching into the preceding block remains a separate route to
13 CNOTs for the complete diagonalizer.
