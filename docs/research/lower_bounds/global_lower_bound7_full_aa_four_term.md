# No real four-term selected image in the adjacent/adjacent span

Consider the ordered six-CNOT schedule
\((01),(23),(12),(01),(23),(03)\). Select an input axis on wire 0
that commutes through its first incident CNOT. Its later active
CNOTs are \((01),(03)\), so after the collective output rotation
used in `global_lower_bound7_degree3_centered_gram.md`, its image
lies in the sound 24-dimensional real Pauli span for the
adjacent/adjacent shape.

**Claim (exact bounded support certificate, ledger run 304).** No
Hermitian involution with exactly four nonzero real Pauli
coefficients in this span satisfies all single-axis spectral-path
identities for the four-cycle shift.

The script enumerates every one of the \(\binom{24}{4}=10,626\)
four-string supports. For \(A=\sum_{r=0}^{3}x_rP_r\), with all
\(x_r\ne0\), a commuting pair contributes a nonidentity Pauli
cross term to \(A^2\). If no disjoint pair produces the same Pauli
string, its coefficient cannot cancel. This immediately excludes
10,176 supports. The remaining 450 supports contain either only
pairwise anticommuting strings or precisely matched commuting-pair
cross terms; the script imposes every resulting cancellation
equation. Unit normalization can be imposed after solving the
homogeneous equations, so it does not change projective feasibility.

On each retained support, restrict the exact integer
spectral-path Gram from run 245 to its ten quadratic monomials.
Full rank modulo the prime 1,000,003 implies full rank over the
rationals and excludes 407 supports. For the other 43, exact rational
rank and Gröbner calculations combine the spectral-path and
involution-cancellation equations, setting \(x_0=1\). This is
lossless because all four coefficients are assumed nonzero.

| Stage | Supports remaining |
|---|---:|
| All four-string supports | 10,626 |
| Cross terms can cancel | 450 |
| Spectral Gram singular modulo 1,000,003 | 43 |
| Joint spectral/involution ideal nonunit over \(\mathbb C\) | 4 |
| Joint ideal has a real nonzero solution | **0** |

The four complex-feasible supports are

```
XXII, XYII, YXII, YYII
XXII, XZII, ZXII, ZZII
XXIZ, XYIZ, YXIZ, YYIZ
YYII, YZII, ZYII, ZZII
```

Each has the same exact projective Gröbner basis,
\(u_1-u_2,\ u_2^2+1,\ u_3+1\). Thus every algebraic point
requires \(u_2=\pm i\), whereas Hermiticity requires real Pauli
coefficients. None is physical. The script also checks saturation
by \(u_1u_2u_3\), so the four surviving ideals are genuinely
nonzero-coefficient complex branches; they are rejected solely by
the real-coefficient condition.

Reproduce with

```bash
python3 archive/legacy_root/global_lower_bound7_full_aa_four_term.py
```

The exact per-survivor ideals, trace coefficient rows, and aggregate
counts are in `global_lower_bound7_full_aa_four_term_result.json`.
The run was bounded to at most 250 modular-singular supports; it
encountered only 43 and completed all 10,626 supports. No real
candidate reaches the cycle trace-moment stage, so those constraints
are vacuous for this branch.

Together with runs 294 and 281, this proves that for **this fixed
representative schedule**, any selected first-axis image of a
diagonalizer must have at least five nonzero Pauli terms in the
fixed collective gauge. It does **not** exclude the whole schedule:
images with five or more terms remain unclassified. The sound
global six-CNOT survivor count remains **3,208**.
