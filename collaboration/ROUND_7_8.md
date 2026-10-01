# Label-free curvature and residual-label escape

Initial tree clean at2ec07b1, three owned workers confirmed completed in the
live registry, no optimizer/certificate process, ledger376 zero running,
board29–32 terminal. Previous goal turn is progress because the improved
invalid12CX lead changed the next experiment. Full goal remains an exact
construction attaining a proved optimum, with real coordinated exchanges;
the accepted13CX upper bound and6CX lower bound do not prove completion.
Details in `work/root_curvature/RECONCILIATION.md`.

## Run377: curvature branch

New board33–35 assigned algebraic interpretation, independent checking and
one numerical experiment. Source is run376 joint0, full12CX+56U3 on14 local
layers with identity optimizer entangler slot1. All168 coordinates remain
free. Budget:100L-BFGS iterations, one168x168 automatic-differentiation
Hessian of label-free ||offdiag(UVU†)||F²/16, then only if smallest eigenvalue
is below-1e-7 two +/-0.1 normalized least-eigenvector starts,150iterations
each;120seconds total including Hessian, one thread, no randomness. Hashes,
schedule and complete code/config were frozen and reserved before execution.

Base fitting converges in47iterations/51evaluations. Actual1.865857092s.
Independent offdiagonal loss0.009423003567534258, maxoff0.1523762940765612;
the loss decreases slightly while maxoff is slightly worse than the source.
It is still invalid and has12CX. Gradient norm2.032059672e-7, Hessian symmetry
error3.33e-16, least eigenvalue-3.343255264e-8 and eigenpair residual4.85e-16.
The threshold is not crossed, so no perturbation fit executes. Hessian
calculation differentiates the floating64 computational graph with inherited
axis-radius regularizer1e-24; it is not an exact arithmetic certificate.

Independent NumPy exponential Pauli-axis matrices, explicit CX bit columns,
and the right-cycle permutation reproduce the saved loss and unitarity. Along
the saved unit eigenvector, central-difference curvature is+4.74846e-8 at
h=1e-3 and-1.34059e-8 at h=5e-4, with objective differences around1e-14.
These checks locate a nearly flat direction; its sign is not robustly resolved
at those finite steps. No certified negative curvature, stationary minimum,
topology exclusion or new lower bound follows. Parameter gauge paths satisfy
f''=vᵀHv+gradient·acceleration; a gauge-mode zero interpretation is clean only
at genuine stationarity. See algebraic and verifier `curvature/REPORT.md`.

## Distinct follow-up authorized after review

The dominant residual pair of the improved lead is(10,14). After reviewing
run377 and the independently checked threshold branch, a separate finite
experiment is reserved to swap these two target labels from the polished base.
Budget: one multiplicity-preserving fixed target,200L-BFGS iterations,
at most five SVD corrections rcond1e-6 with eight backtracking trials,
120seconds total, one thread, no randomness. It is not an additional Hessian
escape or an unchanged restart from the old deletion source.

Algebraic review identifies a potential duplicate: rows11 and15 both carry
label0, so adding swap(11,15) does not change the target. The executed proposal
therefore uses(10,14) alone. Rows10/14 have different labels; swapping preserves
root counts(6,3,4,3). The prior single-wire reduced-target invariant cannot
prove inequivalence to the original run363 base for this new target, and no
such claim is made. Exact certificate work remains separate from numerical
gate-list validation.

## Run378: independently invalid, different metric tradeoff

New board36–38 owns the target derivation, independent review and numerical
experiment. It executes one frozen target in3.655261161s, converging after100
iterations/108evaluations, followed by five accepted SVD corrections, rank60.
Independent12CX list: maxoff0.14301649096287883, offdiagonal loss/16
0.018068152786165515, selected-target loss/16 0.018175553524735204,
maxUV-DU0.08808977407679722. Largest entry improves, but total loss is worse
than the source0.009423003567534258; do not call it a uniformly better basin.
Keep both leads: `work/numerical_label12/curvature/base.json` for lower loss and
`work/numerical_label12/residual_swap/target0.json` for lower maxoff among these
compared saved leads. Neither is a circuit witness or exact upper bound12.

Independent chronological source/target matrix, capacity labels, topology,
source/code/helper hashes and phase reconstruction9.68e-16 agree. Exact
partial trace onto(q1,q3) changes eigenvalues{2,-2,1-i,1+i} to{0,0,1-i,1+i},
proving target inequivalence to the immediate source under product-local output
conjugation. It does not prove inequivalence to the older run363 target,
whose reduced spectrum coincides. No hidden output permutation gate is added.

Both runs are failed/inconclusive, with no exact-certification target or
incumbent promotion. Exact13 and6 <= C_min <= 13 remain; global optimum is
unresolved. Costs are separate from the native set{CX,F(a,b)=exp(i a XX+i b YY),
arbitrary locals} with independently tunable real a,b: preserved6CX+4F,
native depth8, compiled-CX upper14, numerical native JSON plus analytic
regrouping provenance. No direct exact native certificate or10CX claim.

Root board39 derives direction reversal equivalence through local Hadamards.
With arbitrary independent local layers, orientation-only mutations preserve
the realizable family; historical directed hashes are not rewritten and no
canonical family count is asserted. See `work/root_curvature/DIRECTION_EQUIVALENCE.md`.

Recommended next round, not dispatched: one newly frozen architecture run
using the LOWER-LOSS run377 base, two changed unordered-pair schedules rather
than more label swaps. Replace optimizer slot6 pair01 with03 for the first,
or slot8 pair12 with02 for the second; retain identity slot1 and other pairs.
Check historical duplication before reserve; do not substitute direction-only
changes. One warm start per schedule,200iterations each,120seconds total,
one thread, no extra polishing or random seeds. Independent gate-list checking
and exact certification of any passing list remain separate requirements.
These are hypotheses, not promises that the new families admit12CX solutions.

Actual exchanges and actionable board handoffs are preserved in `CHATS.md`
and each researcher directory. All prior immutable evidence remains preserved.
Focused new checks are independent finite differences, source/list/label audits,
script compilation, ledger/board integrity and diff review; shared checker code
is unchanged, so the previously passing22-test infrastructure suite was not
needlessly repeated. Final audit details are in `work/root_curvature/VALIDATION.md`.
