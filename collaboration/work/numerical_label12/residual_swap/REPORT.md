# Single dominant residual label swap

Board38, ledger378, parent377. Source is `curvature/base.json`, with frozen
source/code/helper hashes in config and ledger snapshot. All168local
coordinates remain free on unchanged12CX slot1-I optimizer schedule. Saved
list has56U3+12CX, q0 MSB, chronological right cycle, arbitrary locals and
directed all-to-all CX, ancilla-free. No prior artifact or shared code altered.

Before computation reserved one nearest-root multiplicity-constrained
Hungarian assignment, roots(1,i,-1,-i), counts(6,3,4,3), then ONLY distinct-label
swap(10,14). No(11,15), Hessian reuse, randomness or additional targets.
Budget200L-BFGS iterations/4000evaluations, max5truncated-SVD corrections
cutoff1e-6/eight halving trials each,120seconds overall1Torch/BLASthread.

Runtime3.655261161seconds. Actual100iterations/108evaluations, five accepted
SVD corrections, rank60, alpha0.125. Saved labels
`[2,0,0,3,3,1,2,1,0,2,2,0,1,3,0,0]`. Final gate-list metrics independently
verified:

- maxoff(UVU†)0.14301649096287883;
- unrestricted offdiagonal squarednorm/16:0.018068152786165515;
- selected-target ||UV-DU||F²/16:0.018175553524735204;
- selected-target maximum residual entry:0.08808977407679722.

Candidate invalid. Maximum offdiagonal entry decreases, while Frobenius
objective worsens versus parent3770.009423003567534276. It is not a uniform
improvement; retain distinct best metrics and their sources. Source optimizer
roundtrip agrees upglobalphase9.682678848994767e-16. Read-only diagnostics
retain final serialized-gate metrics; optimizer.loss is before corrections.

Adversarial review checked source/target/counts, chronological fixed12CX
serialization, one-target scope, full168freedom, budgets, derivative/SVD
conventions and objective distinction. No confirmed defect. Independent
verifier confirms full matrix/hash/target/phase; no exact work warranted.
Cooperative wall guard can complete one in-progress evaluation after deadline;
actual computation far below120seconds. All research evidence kept, Python
cache ignored. Coordinator owns commit. No additional search launched.

Next recommendation should use a distinct topology or analytic construction
with a new explicit budget; repeated fitting here establishes no exclusion.
