# Phase conditions for relative-phase selector replacements

Board hypothesis 114. This is a read-only phase audit of the proposed
19-CX factorization. No matrix evaluation, symbolic calculation, candidate
synthesis, or gate-list validation was performed. In particular, no specific
three-CX relative-phase Toffoli implementation is certified here.

## Baseline factorization and budget

Use the Bell decoder prefix (2 CX), the difference router R (2 CX), and the
shared diagonal correction B=CS† on A=(q0,q2) (2 CX). After B, the branch
diagonalizer consists of P (the d1-controlled CX(q0->q2), an exact Toffoli at
6 CX), CH on q0 controlled by d0 (1 CX), and CCH on q2 controlled by
d0=0,d1=1 (6 CX). The exact total is

`2 + 2 + 2 + 6 + 1 + 6 = 19 CX`.

Replacing P and CCH by particular three-CX relative-phase variants would
reduce this count to 13, so one more CX must be removed, for example by
absorbing a residual phase into B. The phrase “three-CX relative-phase
variant” does not specify its left/right phase convention or truth table;
the conditions below must be checked against the actual primitive list.

## When a residual diagonal can pass the final controlled-H gates

Let `D=diag(exp(i phi(d0,d1,a0,a1)))` be a residual computational-basis
phase immediately before the final conditional Hadamards. The diagonal
phase commutes through CH(q0) on the branch d0=1 exactly when

`phi(1,d1,0,a1) = phi(1,d1,1,a1) mod 2*pi`

for both values of d1 and a1. It commutes through CCH(q2) on branch
d0=0,d1=1 exactly when

`phi(0,1,a0,0) = phi(0,1,a0,1) mod 2*pi`

for both a0 values. These conditions are necessary and sufficient for a
diagonal D to commute with the specified H blocks: H mixes the paired target
states, and a diagonal 2x2 phase commutes with H iff its two entries agree.
On d=00 no final H acts, so no equality is required there. A phase depending
only on d always passes because the selectors leave d unchanged.

Any diagonal factor applied after all selector gates is harmless: it simply
changes output eigenvector phases and commutes with the final diagonalized
operator. This freedom does not make a phase before a target-basis rotation
harmless. For example, if an RCCX residual `Delta` is conjugated by target
unitaries V to form a relative-phase CCH, the error is `V Delta V†`. It is
not generally computational diagonal. If Delta varies with the target bit,
the conjugation can mix computational roots, including in control branches
where the exact CCH would be inactive. Target-phase constancy on each pair
mixed by V, or a proven cancellation with the other selector residuals, is
needed before treating the result as an allowed output phase.

For a replacement at P, write the exact relation explicitly, such as
`P_rel = D_out P` or `P_rel = P D_in`. If `D_out` is placed after P, its
pre-P equivalent is `D_in=P† D_out P`, which remains diagonal because P is a
computational permutation. It can be merged into B only if, after all
conjugations and commuting steps, the residual is supported on A and no
selector-dependent cross phases remain. Otherwise it cannot be absorbed into
the fixed two-qubit B without adding another operation.

## One-CX saving condition for B

For an A-only diagonal `D_A=diag(exp(i theta_00),...,exp(i theta_11))`, its
nonlocal phase invariant is

`delta(D_A)=theta_00-theta_01-theta_10+theta_11 (mod 2*pi)`.

The existing `B=CS†` has `delta(B)=-pi/2`. A general two-qubit diagonal is
locally implementable with no CX when its invariant is 0, and is locally
equivalent to a one-CX controlled-Z when its invariant is pi; generic values
use two CX. Therefore merging a residual `D_A` into B saves one CX only if

`delta(D_A B)=pi mod 2*pi`, equivalently `delta(D_A)=3*pi/2 mod 2*pi`.

This is a necessary condition for the proposed one-CX replacement of B's
nonlocal phase, not evidence that either relative-phase selector produces
such a residual. Local phases are free in this accounting. The total residual
must also pass the conditional-H commutation conditions above, or else be
shown to cancel exactly before the H gates.

## Actionable audit for the primitive list

For each proposed relative-phase replacement, the verifier should record its
exact relation to the ideal gate (`D_out G` versus `G D_in`), list the
computational phase function, propagate it through the exact remaining
selectors in their true order, and check the two orbit-constancy conditions.
Then test whether the net phase can be moved adjacent to B and whether its
A-only nonlocal invariant is `3*pi/2` relative to B. If the CCH construction
conjugates an RCCX by a target V, inspect `V Delta V†` rather than assuming
Delta remains diagonal. Any failure requires an explicit cancellation or a
different decomposition; the label “relative-phase Toffoli” alone is not
sufficient. Only after this symbolic phase audit should a full gate list be
checked, and it still needs exact diagonalization certification.
