# Bounded sparse classification of the first AA selected image

For the ordered six-CNOT schedule
\((01),(23),(12),(01),(23),(03)\), the selected wire-0 axis has
active pairs \((01),(03)\). The sound 24-dimensional real Pauli
overapproximation for this adjacent/adjacent tail is defined in
`global_lower_bound7_degree3_centered_gram.md`. The full real
solution set of its spectral-path and involution equations remains
unknown.

Ledger run **290** classifies the finite ansatz

\[
A=\frac{1}{\sqrt m}\sum_{r=1}^{m}s_r P_r,
\qquad m\in\{1,2,3\},\quad s_r\in\{\pm1\},
\]

where the distinct Pauli strings \(P_r\) belong to the 24-dimensional
span and anticommute pairwise. Pairwise anticommutation makes
\(A^2=I\). For at most three distinct Pauli strings it is also
necessary: products from different unordered pairs cannot coincide,
so no commuting cross term could cancel another. We fix the first
sign to +1 because \(A\) and \(-A\)
represent the same unoriented input axis. The exact integer
spectral-path Gram from run 245 tests every such support/sign choice.

| Number of terms | Solutions modulo overall sign |
|---:|---:|
| 1 | 0 |
| 2 | 0 |
| 3 | 4 |

The four solutions are

\[
\begin{aligned}
&(X_0Z_3+Y_0X_1Z_3+Z_0X_1)/\sqrt3,\\
&(X_0Z_3-Y_0X_1Z_3+Z_0X_1)/\sqrt3,\\
&(X_0Y_1Z_3-Y_0Z_3-Z_0Y_1)/\sqrt3,\\
&(X_0Y_1Z_3+Y_0Z_3+Z_0Y_1)/\sqrt3.
\end{aligned}
\]

They are the canonical witness from run 253 and its images, up to
overall sign, under collective \(R_z(\theta)^{\otimes4}\) for
\(\theta=0,\pi,\pi/2,-\pi/2\). These rotations commute with
the cyclic shift. Every solution has
\(|\operatorname{tr}(V_4 A)|^2=4/3\) and
\(\operatorname{tr}(V_4^2 A)=0\).
Thus the sparse ansatz yielded no branch with a different trace
invariant; the previous fixed-A joint obstruction transports to
these four solutions under the same collective rotation and sign.

Reproduce with

```bash
python3 archive/legacy_root/global_lower_bound7_full_aa_sparse_scan.py
```

The exact finite output is
`global_lower_bound7_full_aa_sparse_scan_result.json`. This is a
classification **only** of equal-magnitude one-, two-, and three-term
pairwise-anticommuting sums. It does not classify unequal
coefficients, sums with four or more terms, or the physical two-CNOT
reachable family. The earlier 24-dimensional Gram forced only two
coefficients to zero, leaving 22 real variables before nonlinear
constraints. No six-CNOT schedule is excluded; **3,208** remain.

A concrete next step is to impose the common first-partner Bloch
axis from run 248 on the 22-variable real system, then use the
collective \(Z\)-rotation gauge to fix one coefficient and solve a
low-degree branch decomposition. Each resulting first-axis branch
would need its own second-axis centralizer and joint trace check.
