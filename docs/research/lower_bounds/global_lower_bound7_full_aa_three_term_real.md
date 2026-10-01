# All real selected-axis images with at most three Pauli terms

For the ordered six-CNOT schedule
\((01),(23),(12),(01),(23),(03)\), fix the selected wire-0
axis with active CNOT pairs \((01),(03)\). Its output lies in the
24-dimensional adjacent/adjacent real Pauli overapproximation from
`global_lower_bound7_degree3_centered_gram.md`. Ledger run **294**
classifies every Hermitian involution in this span supported on at
most three distinct Pauli strings that also satisfies the
single-axis spectral-path identities. **Arbitrary real nonzero
coefficients** are allowed.

For at most three distinct Pauli strings, an involution requires all
nonzero terms to anticommute pairwise. Indeed, a commuting pair
contributes a nonidentity Pauli cross term to \(A^2\); products from
different unordered pairs of three distinct strings cannot coincide,
so that cross term cannot cancel. Conversely, for pairwise
anticommuting strings, \(A^2=I\) amounts to unit Euclidean norm of
the coefficients.

For each support, the script restricts the exact integer
spectral-path Gram of run 245 to its quadratic monomials. A full-rank
matrix modulo 1,000,003 is full rank over the rationals and excludes
the support. For each singular matrix it computes an exact rational
rank, fixes the first nonzero coefficient to 1, and solves the
remaining projective quadratic equations with a rational Gröbner
basis. This loses no nonzero real solution because the equations are
homogeneous.

| Terms | Pairwise-anticommuting supports | Full rank modulo 1,000,003 | Singular but infeasible | Feasible supports |
|---:|---:|---:|---:|---:|
| 1 | 24 | 24 | 0 | 0 |
| 2 | 156 | 154 | 2 | 0 |
| 3 | 344 | 330 | 12 | 2 |

The two remaining supports and their complete projective ideals are

| Support in coefficient order | Rational Gröbner basis |
|---|---|
| `XIIZ, YXIZ, ZXII` | \(u_0^2-1,\ u_1-1\) |
| `XYIZ, YIIZ, ZYII` | \(u_0-u_1,\ u_1^2-1\) |

Thus all real solutions on these supports have equal coefficient
magnitudes. Normalization and overall-axis sign leave exactly the
four images listed in `global_lower_bound7_full_aa_sparse_scan.md`.
They are collective \(R_z\) quarter-turn images, up to sign, of
the canonical image A from run 253. The exact joint trace
obstruction for that canonical A, proved in
`global_lower_bound7_joint_aa_family_trace.md`, transports under
these collective rotations because they commute with the shift;
overall sign simply reverses the selected input axis. Therefore
**no first-axis image supported on at most three Pauli strings can
occur in a diagonalizer with this schedule**.

Reproduce the exact support classification with

```bash
python3 archive/legacy_root/global_lower_bound7_full_aa_three_term_real.py
```

The complete surviving support ideals and aggregate counts are in
`global_lower_bound7_full_aa_three_term_real_result.json`. The
excluded full-rank and unit-ideal decisions are reproducible from
the script and the stored run-245 Gram. This theorem covers a
bounded sparse branch of the first-axis 24-dimensional
overapproximation; images with four or more Pauli terms remain
unclassified. It excludes **no entire six-CNOT schedule**. The
sound global screen still has **3,208** six-CNOT survivors.
