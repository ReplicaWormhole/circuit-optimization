# Run402 independent post-fit audit

The initial and best 98-coordinate vectors match their exact seeded/source records and all four saved native/compiled gate-list hashes. Independent chronological matrix checks found coordinate/native/compiled matrices equivalent up to global phase at errors below `8e-16`; unitarity errors are below `1.2e-15`. Both native lists contain 4 CX plus 4 XX/YY gates, and both compiled lists contain 12 CX.

Independent normalized off-diagonal losses are `0.8181143009791337` initially and `0.3750000000000020` at the best point. The max off-diagonal entry moves from `0.9018603550676626` to `0.9999999999999974`; the best circuit is invalid as a diagonalizer. This illustrates that the Frobenius objective improves while the max-entry diagnostic worsens.

The run recorded 127 evaluations, 127 history entries, 113 iterations, and a best point at evaluation 127; the best loss matches both its saved recomputation and minimum history row. It used about 1.365 seconds of the 110-second internal budget, one thread, and at most one final objective recomputation outside maxfun. The fit is bounded numerical evidence only. It produces no exact certificate, incumbent promotion, lower bound, or family exclusion.

The full independent audit is retained in `POSTFIT402_AUDIT.json` here and in `experiments/runs/402/`. No further optimizer, derivative, assignment, or search work was performed.
