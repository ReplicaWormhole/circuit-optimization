# Independent verification of fresh full98 run 393

Hypothesis 76 was claimed by `verifier98fresh`; numerical hypothesis 75 and
ledger run 393 supplied the one-start fit. This work used NumPy/SciPy dense
matrix construction only. It ran no optimizer, derivatives, or perturbation
fits.

## Frozen start and prefit approval

The run's SHA-bound numerical script is
`b19f1aaae55d705fb2487ecb04ca1f8e878688404c2e9a03133edb4e73b305aa`; the
configuration SHA-256 is
`d347c2aeff25828c476354e10baeb0df184022eb8c1ecfc859e7bcc2e7b89b52`.
The 98-vector exactly equals
`default_rng(9261501).normal(0, 0.3, 98)`. Both retained interior rotations
are nonzero: `||A1||=0.43494048`, `||B1||=0.13997853`. Source inspection
confirms all 98 coordinates flow into the family matrix and L-BFGS-B receives
no mask or bounds. The configured limit was one start, one thread, at most
300 iterations and 10,000 evaluations, with 110-second internal and
120-second external limits.

Before the fit, the independent matrix evaluator matched coordinates to the
serialized initial native list to `6.00e-16` up to global phase and the native
list to its compiled list to `4.00e-16`. It counted 4 CX and 4 XX/YY blocks in
the native list and 12 CX in the compiled list. The maximum initial unitarity
error was `8.88e-16`. The frozen chronology was
`L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`.
The prefit details and source hashes are in [PREFIT_APPROVAL.json](PREFIT_APPROVAL.json).

## Independent postfit result

Run 393 stopped after 284 iterations and 301 objective/gradient evaluations
in 3.439 seconds. SciPy reported success because the relative function-change
criterion was met, while the saved best gradient norm `3.90e-8` remained above
the configured `gtol=1e-10`; this is not evidence of stationarity. Its loss
fell from `0.9085013606881812` to `0.4548985707311256`.

I independently rebuilt both complete native and compiled gate lists from
their JSON files for the initial and best vectors. For both vectors, the
coordinate matrix agrees with the native list within `1.1e-15` up to phase,
and the native and compiled matrices agree within `7.5e-16`. Every list is
unitary to at most `1.6e-15`. Each native list has 4 CX plus 4 XX/YY blocks,
and each compiled list has 12 CX. Independent normalized losses agree with
the saved values within `5.6e-16`.

The best compiled list is not a diagonalizer: its maximum target off-diagonal
entry is `0.6743630736` and normalized loss is `0.4548985707`. The max entry
grew from `0.5295689784` even as the Frobenius objective fell. This is one
bounded unsuccessful coordinate-wise fresh start. It excludes neither this
family nor 12-CX circuits, gives no lower bound, and provides no exact
certificate. The established 6-to-13 CX interval and exact13 incumbent are
unchanged.

The full fit payload is in [run 393 result](../../../experiments/runs/393/result.json),
and the complete independent gate audit is in [POSTFIT_AUDIT.json](POSTFIT_AUDIT.json).
