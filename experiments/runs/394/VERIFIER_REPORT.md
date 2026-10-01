# Independent audit of full98 batch run 394

Board hypothesis 78 covers the independent preflight and endpoint audit. I
used the NumPy/SciPy coordinate and gate-list evaluator preserved in
`collaboration/work/verifier98fresh/independent_audit.py`; this batch audit did
not call the optimizer, compute derivatives, or perturb any coordinates.

## Frozen inputs and preflight

The reserved source SHA-256 is
`c8e57295372068ad9ac01535c337d831c49bda0fdc64fcaf5ea31e1285cb7bcd`; the
configuration SHA-256 is
`cc6334176092ff9fd128aa1bb7710a064d297f74f07eb8ef8007204ffb8e3481`. The
dependency hashes matched. All 32 seeds, 9261600 through 9261631, reproduce
the frozen rule: draw 90 local/interior coordinates from a normal distribution
with sigmas cycling through 0.3, 1, 2, 3, then draw 8 F angles uniformly from
`[-pi/2, pi/2]`. The report records a hash and nonzero A1/B1 norms for every
start.

Before the batch, I independently serialized and simulated one initial point
at each scale (seeds 9261600–9261603). Each produced 4 native CX gates and 4
XX/YY blocks, with a 12-CX compilation. Coordinate-to-native phase error was
at most `1.04e-15`, native-to-compiled error at most `7.45e-16`, and native
unitarity error at most `1.56e-15`. The source implements one sequential
single-threaded L-BFGS-B fit per seed, up to 800 iterations and 20,000
evaluations per start, within the 285-second internal and 300-second external
batch caps. The `1e-16` threshold only stops on a promising numerical hit; it
does not certify a solution. Full details are in
[PREFIT_APPROVAL.json](PREFIT_APPROVAL.json).

## Batch and endpoint results

All 32 starts completed in `117.314` seconds; the run did not hit the early
threshold. I independently checked every seed's initial vector against its
RNG rule and preflight hash, then evaluated all 32 saved best vectors against
their complete native and compiled gate lists. I also checked the saved
per-seed hashes, each row's saved objective, and that the global best files
are byte-for-byte copies of the minimum-loss endpoint.

For every endpoint, the native list contains 4 CX and 4 XX/YY blocks; the
compiled list contains 12 CX. Across all 32, coordinate-to-native matrix
error is at most `1.20e-15`, native-to-compiled error at most `6.50e-16`, and
unitarity error at most `2.11e-15`. Independently recomputed normalized losses
agree with the saved values within `2.23e-15`.

The best endpoint is seed `9261616`, with normalized loss
`0.12500000000000003` and maximum target off-diagonal entry
`0.7071067812`. The next best endpoints also settle near loss `0.125` and
maximum off-diagonal entry `1/sqrt(2)`. None is a diagonalizer; none approached
the `1e-16` early-stop threshold. This batch supplies bounded evidence about
these 32 initial points only. It does not exclude the 98-coordinate family,
prove a 12-CX lower bound, or provide an exact certificate. No incumbent
promotion is warranted.

The complete per-seed matrix, residual, cost, and hash audit is in
[POSTFIT_AUDIT.json](POSTFIT_AUDIT.json). The numerical source and all run
outputs remain in [run 394](../../../experiments/runs/394/result.json).
