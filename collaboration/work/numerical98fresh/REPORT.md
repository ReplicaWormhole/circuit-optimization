# Fresh full98 single-start fit

Board hypothesis 75; scientific run 393. This report records one L-BFGS-B fit from the frozen `default_rng(9261501).normal(0, 0.3, size=98)` point in the corrected 98-coordinate family. All 98 coordinates are free; A1 and B1 are explicitly nonzero. The frozen schedule is 4 CX plus 4 XX/YY entanglers, with a 12-CX compiled cost. No restart is authorized.

## Frozen provenance

- Search script SHA-256: `b19f1aaae55d705fb2487ecb04ca1f8e878688404c2e9a03133edb4e73b305aa`
- Configuration SHA-256: `d347c2aeff25828c476354e10baeb0df184022eb8c1ecfc859e7bcc2e7b89b52`
- Reserved run workspace: `experiments/runs/393/`
- The coordinator recovered the already reserved run 393 after a reservation parser error; there was one ledger reservation and one execution.
- Independent verifier reviewed all four initial/best gate lists after execution.

## Fit result

The single-threaded fit ran 284 iterations and 301 objective/gradient evaluations in 3.439 seconds, with no restarts. L-BFGS-B reported success by relative reduction of the objective (`CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH`), not by meeting the requested gradient tolerance: the best gradient norm is `3.8973e-8` versus `gtol=1e-10`. The normalized Frobenius objective fell from `0.908501360688` to `0.454898570731`.

The compiled 12-CX checker reports `valid_diagonalizer=false` and maximum off-diagonal entry `0.674363073634` for the best point. The initial point's maximum off-diagonal entry was `0.529568978351`, so the Frobenius objective improved while this distinct maximum-entry diagnostic worsened. Both gate lists are invalid diagonalizers. The best compiled list is unitary to `7.2e-16`; unitarity does not imply diagonalization.

Independent verifier review found coordinate-to-native phase errors `6.00e-16` initially and `1.07e-15` at the best point; native-to-compiled errors were `4.00e-16` and `7.45e-16`. All four gate lists had unitarity errors at most `1.6e-15`, with counts exactly 4 CX + 4 XX/YY natively and 12 CX compiled. Independently recomputed losses match the saved values within `5.6e-16`. The best list remains invalid, and the ftol stop does not meet `gtol`.

## Interpretation

This is one bounded numerical fit for one fixed 4-CX+4-XX/YY schedule. A coordinate-wise fresh seed does not establish a different gauge orbit or optimizer basin. The improved Frobenius objective and small gradient do not establish a local minimum, an exact diagonalizer, or a global bound. The invalid best gate list does not exclude the family. No exact certificate or incumbent acceptance is claimed.
