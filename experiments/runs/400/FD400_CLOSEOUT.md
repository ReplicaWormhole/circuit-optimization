# Run 400 closeout

The frozen PLAN command was run exactly once and exited 0 in 0.0667 s. It completed all ten independent NumPy coordinate-matrix evaluations under the 30 s, one-thread cap. The saved result is `experiments/runs/400/fd_result.json`; its SHA-256 is `0687bd63036fb928046c43b4cb28349bef779a4fedd8cbdacac58a8a742167f9`. Execution details and all raw samples are retained in `FD400_EXECUTION_AUDIT.json` here and in the run workspace.

The center values reproduce the run398 base objectives: off-diagonal `0.13572330470336316` and fixed-label branch `0.1464466094067264`. Central first differences are close to zero and the saved directional gradient dot products are also at roundoff scale. Central second differences are positive and decrease by about two orders of magnitude when the step decreases by one order, consistent with an O(h^2) finite-step truncation contribution. Finite-step effects and cancellation limit resolution at the saved least eigenvalues near `-1.5e-9`; the tiny AD eigenvalue signs are not independently resolved by these samples. The frozen `-1e-6` escape trigger did not fire. No optimizer fit follows from this run. It supports no local-minimum, global-optimum, topology-exclusion, or exact-certificate claim.

No candidate was generated or promoted. Inputs were not edited; the reconstructed expected run399 audit still matches its historical hash.
