# General opposite-pair obstruction for input-axis images

Let V be the four-qubit right shift and U a diagonalizer with U†VU=D
computationally diagonal. Let A=U p_j U† for any unit input Pauli axis.
Then A is a traceless Hermitian involution and satisfies all distinct-label
spectral-path identities P_a A P_b A P_c=0.

Suppose A has support on physical opposite wires (0,2). Expand its
traceless two-qubit factor H in the 15 nonidentity Pauli products with
real coefficients x_ab. The stored exact spectral-path Gram G has
rational rank75 and nullity45. For each of the six unordered distinct
Pauli-factor pairs, the linear monomial functional

    x_ab² + x_ba² - 2 x_ab x_ba

annihilates every exact rational nullvector of G. The certificate script
checks all270 exact rational equalities. Spectral-path zero puts the
monomial vector in ker G (because G is a positive semidefinite Gram).
Hence each displayed expression is zero, and reality implies x_ab=x_ba.
Therefore SWAP02 H SWAP02=H. In particular A commutes with V²,
which equals SWAP02 SWAP13.

The two-qubit swap has a three-dimensional symmetric space and a
one-dimensional antisymmetric singlet. Since H commutes with swap, its
singlet eigenvalue h is real; H²=I implies h=±1. Tracelessness gives
tr(H|symmetric)=-h. Consequently

    tr(SWAP02 H)=-2h,   tr(V² A)=2 tr(SWAP02 H)=-4h ≠ 0.

However, any proper-support input-axis image commuting with V² has
tr(V² A)=0: if p_j=±Z_j then A commutes with V and must have full
support, a contradiction. Otherwise p_j has a nonzero off-diagonal
entry. From [p_j,D²]=0 the D² entries agree on computational strings
paired by bit j; their contributions to tr(D²p_j) cancel. Cyclicity
of trace then yields zero. This contradicts -4h.

Thus no input-axis image of any diagonalizer can have support confined
to an opposite physical pair. Cyclic relabeling gives pair(1,3) too.
Combined with the earlier arbitrary-adjacent-pair obstruction, any
proper input-axis image confined to two wires is excluded.

# Exhaustive necessary schedule screen

For each input wire, choose the Pauli axis that passes unchanged
through its first incident CNOT (control Z or target X, pulled back
through preceding local rotations). Its post-gate support is one wire.
Subsequent local gates do not enlarge support; CNOT directions do not
affect the support-cone overapproximation. The new screen excludes
any prior survivor having a two-wire opposite cone for this selected
axis. All46656 six-CNOT unordered-pair schedules are enumerated.

Of the prior3208 survivors,184 are newly excluded and3024 remain:
1200 of degree3333,1680 of degree4332,144 of degree4422. All40
previous degree5322 survivors are eliminated. The first excluded
representative is (01),(02),(23),(02),(13),(13), with wire1 and wire3
selected cones confined to opposite pair(13).

This is a necessary screen, not a seven-CNOT lower bound. The interval
remains6≤Cmin≤14. Surviving schedules are not circuit witnesses.

Reproduce:

    python3 global_lower_bound7_opposite_pair_support_certificate.py
    python3 global_lower_bound7_opposite_pair_support_audit.py

Ledger runs356 and357 archive both scripts immutably. The independent
audit reconstructs the entire Gram using Gaussian integer int64
arithmetic and matches every entry. Block coordinates are at most16;
Gram sums have12288 terms and bound3145728, so no int64 overflow occurs.
A separate bitmask cone implementation agrees for every enumerated
six-CNOT schedule. Root independently checked the rational nullspace
and all six squared-difference functionals.

Adversarial review: no confirmed issue. The proof uses real coefficients
essentially, justified by Hermiticity, and the trace argument uses the
involution assumption essentially, justified by input Pauli conjugation.
No assertion is made about arbitrary non-involutive symmetric operators.
Existing files and user changes were preserved; no commit was made.
