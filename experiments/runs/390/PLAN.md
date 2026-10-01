# Run 390: independent full98 directional finite difference

This is one independent SciPy/NumPy matrix objective audit along the least-eigenvalue direction from full98 Hessian run389. It was separately reserved in the scientific ledger before any objective evaluation and linked to board hypothesis 70. Run389 reports lambda_min=-2.6985205738843378e-9, which fails the -1e-6 escape trigger; the coordinator requested this independent derivative check regardless of the trigger. No escape fit or optimizer is authorized.

Frozen input: `FROZEN_INPUT.json`, SHA-256 `92287d09ce0424c0a450cbe7198dbc3aa7027dbc8f4e435e6ee43c55f46d899d`. It contains seed 9261101, the 98-coordinate point, the normalized run389 least eigenvector, steps h=1e-3 and 5e-4, exactly five objective evaluations, one thread, and a 120-second external timeout. The ledger reservation records source run389 result SHA-256 `bd095f709735a600de28b420420bbab1526f2b14bdb64123bd66d088d0988210`, full input SHA-256, source Hessian script hash, code hash, steps, evaluation count, no optimizer, thread count, and timeout. The point and direction hashes are present in the hash-bound FROZEN_INPUT file.

Method: independently form the full 16x16 native family matrix using NumPy/SciPy rotations, CNOT permutations, and matrix exponentials for XX+YY gates. Evaluate normalized cycle off-diagonal loss at x and x +/- h v, then use the centered second difference. The external invocation uses `timeout --signal=TERM --kill-after=5s 120s` with OMP, OpenBLAS, MKL, and NumExpr limited to one thread.

Interpretation is limited to this one directional finite-difference check of the AD Hessian. It does not certify a full Hessian, local/global minimum, diagonalizer, topology exclusion, or exact circuit.
