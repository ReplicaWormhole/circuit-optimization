# Dynamic eigenlabel full98 fit review

Board hypothesis 83. This is a read-only source/config review of the coordinator's next single-start full98 eigenlabel fit. No objective or matrix evaluation, optimizer, derivative, or certificate is run here.

## Review checklist

- Confirm the assignment cost matrix matches `min_D ||U V4 - D U||_F^2/16` with the exact repeated spectrum `+1` x6, `-1` x4, `+i` x3, `-i` x3.
- Confirm the Hungarian assignment uses the frozen deterministic tie behavior and records assignment changes.
- Confirm each gradient is computed with the selected diagonal `D` fixed for that evaluation, and the optimizer is allowed to cross assignment boundaries by reassignment on its next evaluation.
- Confirm the one-start seed/source, all98 free coordinates, limits, and no-restart guards match the ledger reservation.
- Record any code-path or logging issue; do not run any fit or objective check.

Frozen files: `experiments/runs/395/fit.py` and `experiments/runs/395/config.json`.
Script SHA-256 `27a8ff0747dbf79f555fe4c4dec429cde1cc5530644e9b9715f6f2289e4ffa50`;
config SHA-256 `7e1a8e0528f347b3e05e76d4cc9c617f264c0be7d0d9d1273b2c48530f7ac8c2`.

## Static review findings

The label slots have the exact multiplicities `+1` x6, `-1` x4, `+i` x3,
`-i` x3. For each row, the script's raw Hungarian cost is
`sum_s |(UV)[r,s] - lambda U[r,s]|^2`; summing selected row costs is exactly
the full `||UV-DU||_F^2`, including the assignment-independent off-diagonal
part. Its selected assignment agrees with the reduced diagonal cost
`|B_rr-lambda|^2`, but that reduced cost alone would omit the original
off-diagonal loss. The script correctly compares its autograd value with the
full assigned row-cost sum; **do not add the off-diagonal loss again**.

The selected labels are converted to a constant complex Torch diagonal, so the
gradient is the derivative of the selected fixed-D branch. Assignment is
recomputed on every objective call. The history stores the selected slot
columns, both the label and original losses, gradient norms, and evaluation
numbers; the frozen slot order makes labels reconstructible. Permuting slots
that carry the same eigenvalue can change `columns` without changing `D`; only
changes in the row-to-root label vector are substantive branch switches.
Distinct-root ties remain kink points, so L-BFGS-B's smooth-model assumptions
can be unreliable near a switch; the history is sufficient to diagnose this
afterward.

The initial 98-vector is loaded exactly from run394 seed9261620's saved best
point; no random draw or perturbation occurs. All coordinates remain free. The
fit has one start, no restart path, a pre-evaluation 20,000-call guard, the
frozen 1,000-iteration cap, a 110-second monotonic fit deadline, and a
120-second external cap. It stores initial, best-label, best-original,
last-evaluated, and optimizer-return point lists and checker summaries. A
nearhit is recorded before early exit, not called a solution. Script/config
hashes match the frozen values, and an AST parse passed. No objective, matrix,
gradient, optimizer, or certificate calculation was run in this review.
