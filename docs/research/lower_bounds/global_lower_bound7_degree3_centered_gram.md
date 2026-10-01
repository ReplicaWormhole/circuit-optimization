# Bounded 24-dimensional degree-three selected-axis calculation

This is an exact necessary-condition calculation for the six-CNOT
four-qubit cyclic-shift problem. It does **not** exclude a new schedule
or prove a seven-CNOT lower bound.

Suppose an input wire \(j\) meets three CNOTs. Choose an input Pauli
axis that commutes through its first incident CNOT. If precisely two
later CNOTs are active in its potential support cone, they are the
remaining two CNOTs incident to \(j\), on ordered pairs \((j,k)\)
then \((j,\ell)\). Any intervening CNOT disjoint from the cone acts
trivially on this selected image. Up to dihedral symmetry of the
physical four-cycle, with chronological order retained, the distinct
physical triples are:

| Shape | Representative \((j,k,\ell)\) |
|---|---|
| Both pairs adjacent to \(j\) | \((0,1,3)\) |
| Adjacent then opposite | \((0,1,2)\) |
| Opposite then adjacent | \((0,2,1)\) |

After the first active CNOT, the selected image has a nonidentity
Pauli factor on \(j\) and optional factors on \(k\). The second
active CNOT, joining \(j,\ell\), preserves a nonidentity factor on
\(j\); all factors it creates on \(\ell\) share one Bloch axis.
A collective single-qubit rotation commutes with the shift and can
put that axis on \(Z_\ell\). Hence every selected image lies in the
overapproximate 24-dimensional real span

\[
\operatorname{span}_{\mathbb R}
\{\sigma_\mu^{(j)},\;
  \sigma_\mu^{(j)}\sigma_\nu^{(k)},\;
  \sigma_\mu^{(j)}Z_\ell,\;
  \sigma_\mu^{(j)}\sigma_\nu^{(k)}Z_\ell:
  \mu,\nu\in\{X,Y,Z\}\}.
\]

This span allows arbitrary one-qubit gates and either CNOT direction.
It discards further constraints of the actual two-CNOT image, such
as the common axis on \(k\), so failure of the span test would be a
sound obstruction.

`global_lower_bound7_degree3_centered_gram.py` constructs an exact
integer spectral-path Gram for each shape on the 300 quadratic
monomials in these 24 coefficients (ledger run 245). Every entry is
an exactly represented Gaussian-integer result, with maximum
absolute value 81,920. Exact rational nullspace calculations give:

| Ordered shape | Exact Gram rank | Nullity | Coefficient squares directly forced to zero |
|---|---:|---:|---|
| Adjacent, adjacent | 255 | 45 | \(Z_jX_k\), \(Z_jY_k\) with \(Z_\ell\) present |
| Adjacent, opposite | 260 | 40 | \(Z_jX_k\), \(Z_jY_k\) with \(Z_\ell\) present |
| Opposite, adjacent | 236 | 64 | None |

The exact stored basis labels for the forced coefficients are `ZXIZ`
and `ZYIZ` in the first representative, and `ZXZI` and `ZYZI` in
the second. These two forced squares leave 22 real coefficient
variables in each of the first two shapes; the third leaves all 24.
The remaining rank-one equations, combined with the Hermitian
involution equations, were **not** solved in this bounded run. Rank
deficiency alone does not establish a physical image. The latter
two shapes already contain same-axis opposite-pair Pauli products
that satisfy the spectral-path identity, so spectral-path positivity
alone cannot exclude those entire spans. The prior half-shift trace
obstruction excludes each such product as the image of an input
Pauli axis, but a classification of all remaining solutions is
needed before excluding a schedule.

No new six-CNOT schedule was excluded. The earlier sound filter
still leaves **3,208** ordered unoriented pair schedules, including
**1,272** of degree 3333. A survivor is not a verified diagonalizer;
the global interval remains \(6\le C_{\min}\le14\).
