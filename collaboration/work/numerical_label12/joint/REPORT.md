# Two composed label targets

Board31, ledger376, source and12CX topology unchanged from run375. No previous
artifact overwritten. All168 local axis coordinates remain free; identity
slot1 is optimizer-only. All saved candidates have12CX+56U3, chronological,
q0 most significant, ancilla-free directed all-to-all CX plus arbitrary locals.

Frozen before execution: verified base label vector
`[2,0,2,1,3,1,0,3,0,2,0,0,1,3,2,0]`, roots(1,i,-1,-i), counts(6,3,4,3).
Target0 composes(2,6)(3,7); target1 composes(8,13)(9,12). No numerical target
screening or randomness. Independently warm-start each at the SAME saved
run363 slot1 deletion source, SHA9a9dd0e0b430e935f99f4d998fb25d5e1b48160ad16a6c46a27e2574c59c2398.
Budget: at most200L-BFGS iterations/4000evaluations each, five truncated-SVD
Newton corrections each, cutoff1e-6, eight halving trials per correction,
120seconds overall and one Torch/BLAS thread. Config+script snapshot in ledger.

Runtime8.639925215975381seconds. Target0 converged99iterations/111evaluations,
five SVD corrections retained rank60; target1 converged151/167, five
corrections rank59. All corrections accepted. Final serialized-gate checks:

|Target|fixed ||UV-DU||F²/16|max offdiag(UVU†)|valid|
|---|---:|---:|---|
|0|0.009459905009667346|0.1517277915509333|false|
|1|0.46343390751450736|0.7071067918307067|false|

First target improves this saved12CX numerical basin compared with run375
single swaps and source maxoff0.258819; it is still invalid. No reduction,
exact candidate, topology exclusion or global lower bound follows. Fixed
target loss differs from offdiagonal loss. `optimizer.loss` is before Newton;
`diagnostics.json` supplies final full serialized-gate values. Source roundtrip
matches upglobalphase1.4488835837920476e-15. Assigned verifier independently
checks matrices and target provenance; no exact13 recertification.

Independent verifier agrees: unrestricted offdiagonal loss0.00942310839202986
for target0 and0.3750000000000001 for target1. Read-only `residual_pairs.json`
records the exact saved labels and residual entries. Target0 dominant pair is
(10,14), two-direction squared energy0.04604264545784692; then(4,8),(5,9),
(0,12),(1,13), each about0.0133820547. These are diagnostics of an invalid
candidate, not exact constraints.

Adversarial review checked budget, source hash, unchanged chronological12CX
schedule, exact prescribed label permutations, full168freedom, independent
source resets, rejected-step handling and numerical scope. No confirmed issue.
Cooperative wall guard permits one already-started evaluation/Jacobian/SVD to
complete after deadline; actual run is far below cap. Kept all research evidence,
ignored Python cache. Shared checker unchanged; coordinator owns commit.

Recommendation: newly budgeted label-free refinement from this improved
target0 (different initialization from run363), or its largest residual
(10,14) spectral-label swap. Finite perturbation starts require a frozen seed,
amplitude and subspace rule; Jacobian nullspaces may largely describe gauge
freedom. No additional search launched.
