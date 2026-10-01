# Independent verification plan for run396

Board hypothesis 89 is the claimed verifier assignment. A first proposal, board88, depended on still-running executor board87 and could not be claimed under board dependency policy; board89 is the replacement verifier assignment without a dependency. Board87 remains the only executor/run-linked hypothesis. No board or ledger rows were directly edited.

## Preflight

Run396 source, config, and imported assignment hashes matched the reservation: batch SHA256 `c0e5eeed37086b68f588d291c9fdfed69a71043e8fb52646c264bbc7b24c0c47`, config SHA256 `775abe1ab572124817e907db842df1c938a5f4dbb0f2998c0c1aa34fb3a9194c`, and run395 assignment SHA256 `27a8ff0747dbf79f555fe4c4dec429cde1cc5530644e9b9715f6f2289e4ffa50`. SciPy is 1.16.3. The eight starts reproduce exactly from `default_rng(seed)`, with seeds 9261700–9261707, 90 normal coordinates at cyclic sigma 0.3/1/2/3, followed by eight uniform angles on `[-pi/2,pi/2]`. Their full, normal, and uniform vector digests are recorded in `PREFIT_APPROVAL.json`.

Independent NumPy coordinate matrices and native/compiled gate-list matrices were checked for one start at each of the four scales. All have 4 native CX + 4 native XX/YY and 12 compiled CX; coordinate/native and native/compiled phase errors are below `9.5e-16`, and unitarity errors below `1.45e-15`. On all four representatives, independent raw Hungarian costs agree with direct `||UV-DU||_F^2/16` and its conjugated form to within `2.3e-16`; selected multiplicities are 6,4,3,3. No optimizer or derivative was run by this verifier.

The static reviewer notes one additional original-objective diagnostic call after each completed fit is not included in the per-fit `calls` counter for `maxfun`; this does not alter the L-BFGS-B evaluation cap. If the batch deadline expires before a new start, the current stop path is `no_evaluated_point`, while earlier completed rows are saved. Root was notified to include these bounded implementation limitations in the run outcome.

The hash-bound preflight is at `PREFIT_APPROVAL.json` and its exact copy is in `experiments/runs/396/PREFIT_APPROVAL.json`. Post-fit endpoint audit will be appended here after the reserved batch completes.

## Post-fit audit

Run396 completed all eight starts within 42.77 seconds, below the 110-second internal / 120-second external batch cap. The independent auditor checked every row's seed, sigma, exact initial vector, saved gate-list hashes, coordinate-to-native and native-to-compiled matrices, unitarity, cost counts, recomputed Hungarian assignment, label loss, and original off-diagonal loss. It also checked all saved evaluation-history assignment columns have the required 6/4/3/3 slot multiplicities and that each row's best assignment agrees by label vector with the saved best point (duplicate slots can permute columns). The exported global-best lists exactly match the lowest-loss row. All eight are invalid by the native checker; no endpoint is an exact diagonalizer.

Independent best is seed9261703 (sigma3), label loss `0.1464466094067264`, original off-diagonal loss `0.13572330470336316`, and max off-diagonal entry `0.35355339341816655`. This is the same plateau, within roundoff, as run395's seed9261620 warm start. Other label losses by seed are: 9261700 `0.7939717829`, 9261701 `0.5329642289`, 9261702 `0.7037866308`, 9261704 `0.4946323502`, 9261705 `0.7685818400`, 9261706 `0.5287263146`, 9261707 `0.6735567184`. Every serialized list pair passed matrix, unitarity, and cost checks; all saved checker flags are false. The largest independent coordinate/native phase error over all rows is `1.01e-15`; the largest native/compiled phase error is `6.78e-16`; unitarity errors remain below `1.56e-15`.

Result hash: `bf7a4d1d32f3f26c53f353f729a2e35f72a6c636ac106800cb709c7087c286a4`. Full per-seed endpoint measurements are in `POSTFIT_AUDIT.json`, copied to the reserved run directory. This bounded batch provides no exclusion, global-optimality statement, exact certificate, or promotion. The repeated plateau is evidence only about these eight starts and this frozen schedule/objective/budget.
