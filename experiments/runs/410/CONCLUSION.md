# Run410 independent endpoint audit

Passed all twelve endpoint cases and84 full matrix builds in0.1797seconds,
within60s one-thread bound. Replayed exactly1728 normal PCG64 draws for seed
10012101 and verified all144-coordinate base/initial/best vectors against
independent SciPy-expm and primitive-gate constructions, eleven-CNOT counts,
full gate chronology, histories and saved metrics. Maxphase1.78e-15,
maxmetricdiff6.09e-15, maxlossdiff2.00e-15. Ninety-three input artifact hashes
are frozen in input_manifest.json. No optimization or derivative/exact check.

Best9 and second7 remain invalid at1e-9 despite small numerical residuals.
The audit establishes consistency and honest near-hit provenance; no topology
exclusion, local minimum, exact diagonalizer or incumbent promotion follows.
Run closed complete/numerical; all inputs and outputs retained. The subsequent
analytic11 proposal is a distinct by-hand construction needing its own frozen
complete-list exact certification.
