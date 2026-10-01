# Dominant residual label swap: structural scope

Hypothesis36; structural math only, no numerical or symbolic screening.
Source is run377's invalid twelve-CX `numerical_label12/curvature/base.json`.
Numerical researcher owns one bounded fixed-target fit: nearest multiplicity
constrained labels, swap10/14,200iterations plus at most five SVD corrections,
120seconds and one thread; all168 local coordinates free. Roots are i^label,
q0 is MSB, chronological gates and ancilla-free directed all-to-all CX plus
arbitrary local gates define the cost model.

The saved diagonal entries of377 are closest to the same labels as the first
joint target A: [2,0,0,3,3,1,2,1,0,2,0,0,1,3,2,0]. These already have
cycle capacities(6,3,4,3); numerical researcher must nevertheless freeze the
actual constrained assignment from its stated rule. Swapping10/14 exchanges
1 and-1 while preserving capacities. Rows11/15 both have label0, so adding
that second swap changes no target and must not count as a distinct proposal.

A new reduced-target invariant distinguishes this target from A under free
product-local output conjugation. Trace D over q0,q2, retain q1,q3. Its
diagonal sums in order00,01,10,11 are [2,1-i,-2,1+i] for A. Swapping rows10/14
subtracts2 from the first sum and adds2 to the third, giving [0,1-i,0,1+i].
The resulting eigenvalue multisets differ. For product W, partial trace of
W D W† is conjugated by W1 tensor W3, so its spectrum is invariant. Thus
the swapped target is not product-local equivalent to A. The historical
pre-joint base has the same reduced spectrum as the swapped target; neither
this nor the earlier one-wire invariant proves inequivalence to that base.
Other gate-cost-preserving symmetries remain unclassified.

This selects a fixed spectral objective, not a hidden appended permutation
gate. Every diagonalizer may use arbitrary output-label order and arbitrary
orthonormal bases inside root eigenspaces. Failed fitting at one fixed target
does not exclude this topology or other targets/bases. Report selected-D loss,
label-free offdiagonal loss, chronological count, unitarity and saved hash.
Any passing candidate needs independently validated gates and independent exact
certification before incumbent promotion. The exact13CX incumbent persists.

Run377's least Hessian value approximately-3.34e-8 did not satisfy the frozen
-1e-7 trigger. Independent finite differences had unresolved nearzero sign,
and gradient norm approximately2e-7. No curvature escape was authorized by
that rule, and neither local minimum nor saddle was certified. The label-swap
round is a distinct mechanism, with its own reservation and budget.

Adversarial review: multiplicities preserved; reduced traces conserve total
trace2; no local permutation treated as free merely because output labels
are unconstrained mathematically. No confirmed defect. Remaining risk is
heuristic fit success and incomplete symmetry classification. Files are
intentional kept artifacts; coordinator owns the shared commit.

Verifier independently reproduced both retained two-wire diagonal sums using
direct Gaussian-integer arithmetic and confirmed the stated invariant scope.
