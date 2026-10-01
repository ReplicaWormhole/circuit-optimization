# A balanced-cut Schmidt-rank obstruction

This is a universal necessary condition for an ancilla-free
diagonalizer \(U\) of the four-qubit cyclic shift \(V_4\), allowing
arbitrary one-qubit gates, CNOT directions, eigenvalue orderings, and
orthonormal bases inside degenerate eigenspaces. It does **not** prove
the global seven-CNOT lower bound.

## Exact eigenspace calculation

Let \(E_i=\ker(V_4-iI)\), a three-dimensional eigenspace. The three
columns chosen by the exact projector calculation in
`global_lower_bound7_cut_schmidt.py` have disjoint computational-basis
support, squared norm four, and span \(E_i\). Divide them by two to
obtain an orthonormal basis \((b_0,b_1,b_2)\). Write any vector in
this space as \(\psi=x_0b_0+x_1b_1+x_2b_2\).

Across the physical cut \(01|23\), reshape \(\psi\) into a \(4\times4\)
coefficient matrix. All its \(3\times3\) minors vanish precisely when

\[
q(x)=x_1^2-2x_0x_2=0.
\]

The same equation holds across \(03|12\), using the same basis. The
script computes the minors with exact Gaussian-rational SymPy
arithmetic. Its Gröbner basis for their ideal, up to nonzero scalar
factors, is \((x_0q,x_1q,x_2q)\). Since \(\psi\ne0\), this proves
\(\operatorname{SchmidtRank}(\psi)\le2\) if and only if \(q(x)=0\).
Across \(02|13\), every vector in \(E_i\) has rank at most two, so
this particular argument gives no crossing bound there.

Put

\[
Q=\begin{pmatrix}0&0&-1\\0&1&0\\-1&0&0\end{pmatrix},
\qquad q(x)=x^{\mathsf T}Qx.
\]

The matrix \(Q\) is symmetric and unitary. Suppose an orthonormal
basis of \(E_i\) consisted entirely of Schmidt-rank-at-most-two
vectors across either adjacent cut. Its coefficient columns would
form a unitary matrix \(R\), so \(S=R^{\mathsf T}QR\) would be a
symmetric unitary \(3\times3\) matrix with zero diagonal. No such
matrix exists: writing its upper off-diagonal entries as \(a,b,c\),
orthogonality of its three pairs of rows gives
\(b\bar c=a\bar c=a\bar b=0\). At most one entry can be nonzero,
leaving a zero row. This contradicts unitarity.

Therefore **every orthonormal eigenbasis of \(V_4\) contains a vector
of Schmidt rank at least three across each of the cuts \(01|23\)
and \(03|12\)**. The vectors witnessing the two cuts need not be the
same.

## Circuit consequence and six-CNOT count

Every column \(U|z\rangle\) of a diagonalizer is an eigenvector of
\(V_4\), and its three columns in \(E_i\) form an orthonormal basis
of that sector. Starting from the product state \(|z\rangle\), a CNOT
crossing a fixed bipartition at most doubles Schmidt rank; a CNOT
within either side, and any one-qubit gate, leaves that rank unchanged.
With at most one crossing CNOT, every column would have rank at most
two, contrary to the preceding theorem. Hence **at least two CNOTs
must cross each adjacent balanced physical cut**.

`global_lower_bound7_cut_crossing_filter.py` applies this exact
condition after the necessary conditions of
`global_lower_bound7_schedule_filter.py` (ledger runs 229 and 231).
All **3,576** six-CNOT schedules left by the earlier screen already
have at least two CNOTs across both cuts. Thus this new theorem
excludes **zero** additional schedules. It provides no improvement
to the established global interval \(6\le C_{\min}\le14\).
