# Run 391: corrected independent full98 directional finite difference

Run390 is preserved as a failed pre-evaluation guard check with zero objective evaluations. Run391 is a distinct board hypothesis (71) and newly reserved ledger row (391), made after the source layout error was diagnosed. It uses the same frozen point and least-eigenvector direction from run389 without renormalization or coordinate change.

The complete frozen input is `FROZEN_INPUT.json`, SHA-256 `57f86dfdb0cf8e7b7d585cc370cc574b8506fdd963adbeb164d5cbb93e72f915`. It binds run389 result hash `bd095f709735a600de28b420420bbab1526f2b14bdb64123bd66d088d0988210`, Hessian script hash `739c8585bcb5b0146a3d735dd06dc2c0be65a7eaf9c923ffffa766303e741957`, exact point/direction JSON hashes, step sizes `1e-3` and `5e-4`, five objective evaluations, one thread, 120-second external timeout, and no fit authorization. The corrected finite-difference script SHA-256 is `9b0c7dd4b1653700ea1cfe1863babb2ea68cf24f5c81591aca7acea611d7dbee`.

Method: independent NumPy/SciPy dense matrices for all one-qubit rotations, CNOTs, and XX+YY gates; evaluate `f(x)` and `f(x +/- h v)` at both frozen steps and form centered curvatures. The saved eigenvector is serialized as a list of columns, so the script checks the first stored column directly. A pre-evaluation check confirmed exact point/direction equality with run389 and unit norm within `2e-12`.

The invocation uses an external `timeout --signal=TERM --kill-after=5s 120s` and OMP, OpenBLAS, MKL, and NumExpr thread counts of one. No optimizer, starts, or escape fits are part of this run. Results only check one Hessian direction and do not prove a minimum, diagonalizer, or topology exclusion.
