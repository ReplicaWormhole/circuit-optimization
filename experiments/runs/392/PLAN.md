# Run 392: independent FD along exact run389 least eigenvector

Run392 is the coordinator-approved corrected run after run390/391 exposed an eigenvector JSON orientation mistake. The saved `eigenvectors_columns` field is a nested 98x98 row-major JSON array; its first column, `Q[:,0]`, is the least eigenvector. Run391 mistakenly used row `Q[0,:]`; those five evaluations are preserved but do not validate the Hessian minimum. Run390 performed zero objective evaluations.

This run uses the exact column `Q[:,0]` with no renormalization, and the exact run389 point. Root independently confirmed array equality, point equality, source/result/code hashes, vector norm `1.0000000000000002`, and least eigenpair residual `6.9848e-16` before reservation. The frozen input is `FROZEN_INPUT.json`, SHA-256 `3e7cbe2460b990942d4463f053fd3377973cf22bc3855b03da64568f3a941118`; code SHA-256 `9102abd85a61e062682361f455c3e27a183f3bf94cb6bda0c3b9d60921e86741`.

The scientific ledger reserved run392 before evaluation and records the run389 result hash, frozen-input hash, point and direction hashes, source Hessian code hash, method, centered step sizes `1e-3` and `5e-4`, five objective calls, no optimizer/iterations/fits, one thread, and an external 120-second timeout. The computation independently forms the 16x16 family matrices with NumPy/SciPy and evaluates `f(x)` plus `f(x +/- h v)` at both steps.

The invocation uses an external `timeout --signal=TERM --kill-after=5s 120s` and OMP, OpenBLAS, MKL, and NumExpr thread counts of one. No other scientific stage is authorized. Results check only one direction and do not establish full local minimality, a diagonalizer, topology exclusion, or an exact certificate.
