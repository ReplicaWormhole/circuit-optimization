# Independent curvature diagnostic verification

## Preflight and run397

The seed9261703 best point from run396 was independently bound to result SHA256 `bf7a4d1d32f3f26c53f353f729a2e35f72a6c636ac106800cb709c7087c286a4` and f64le point digest `d98c9c14ca11d83a31f6682a5b742c9ab45ee2bbe5041ec9d0d98934cc7c1c49`. The saved native/compiled lists match the independent coordinate matrix up to phase at `1.01e-15` and each other at `4.91e-16`; unitarity errors are below `7.8e-16`. Counts are four native CX, four native XX/YY, and twelve compiled CX. Independent finite losses are labeled `0.1464466094067264`, raw assigned cost `0.14644660940672655`, conjugated labeled residual `0.1464466094067265`, original off-diagonal `0.13572330470336316`, and maximum off-diagonal entry `0.35355339341816655`. The candidate is invalid.

Run397’s script/config hashes matched their frozen values, and static review found its branch definitions and row-root exclusions correct. The timeout path could stop in the initial assignment-gap loop before serializing any output. Root therefore closed run397 before execution; no Hessian or gap search was run there. Run398 repaired the timeout handling by putting assignment and Hessian-row loops inside one catch and retaining partial results.

## Run398 preflight and post-fit audit

Run398 code SHA256 is `5e67087e57f9b7d3a8915ac41339db7877b22df958a1c60e67896456cac9b086`; config SHA256 is `a5e6ecbf165e548321c5b178edb0c4d2941cfa72dfbce2179fd89b95acaeae13`. All frozen dependencies and the exact run396 point hash matched. Static checks confirm two objectives in the same 98-coordinate chart, fixed labels from the frozen assignment, 16 row bans excluding every duplicate slot for the selected root plus the base assignment, at most two first gradients and 196 Hessian rows, one thread, and 110/120 second internal/external caps. No AD/full-Hessian or gap calculation was performed for preflight.

The completed output is bound to result SHA256 `8ff1a48ab94146837f3756eb86e36bad87fc9dd2a0bbc78748371eddb651d1c2`. The base gate lists and objective values independently match. I checked all sixteen recorded alternative assignments for row-ban feasibility and recomputed their costs from the independent raw cost matrix without performing more Hungarian solves. The stored gap `0.2133883472436626` agrees with the minimum saved alternative cost minus base cost (`0.2133883472436623`); the auditor does not independently re-solve each constrained assignment.

Both stored Hessians pass independent eigensystem checks. Maximum symmetry error is `4.44e-16`; recomputed eigenvalues differ by at most `6.7e-15`; eigenvector columns are orthonormal to `2.6e-15`; maximum eigensystem residual is `7.3e-15`; all signs follow the frozen largest-entry-positive convention. Least eigenvalues are `-1.4609760402e-9` (off-diagonal) and `-1.5584236709e-9` (fixed-label branch), with stored gradient norms `1.7761e-8` and `1.9164e-8`. Neither reaches the frozen `-1e-6` trigger. This is floating-point curvature evidence at one point, not a local-minimum proof or an exclusion.

The full output audit is in `POSTFIT398_AUDIT.json`, copied into `experiments/runs/398/`. Board96 is closed as the verifier assignment for run398; board95 was accidentally claimed during concurrent ID creation, no work was performed there, and it was closed inconclusive with that reason.

## Finite-difference validation prepared, not run

`fd_check.py` and `fd_manifest.json` freeze exactly the run398 point, fixed labels, and the least-eigenvector column 0 for each Hessian, including full vector hashes. They use steps `1e-2` and `1e-3` and ten independent matrix evaluations: two separate center evaluations and eight signed perturbation evaluations. Each perturbed point is used only for its matching objective; no assignment is recomputed at perturbed points. The script records central first and second differences against the saved gradient projection and least eigenvalue. It has only been syntax-checked. Do not execute until root reserves and board-links run399.

No optimizer, AD/full-Hessian computation, finite-difference evaluation, or extra Hungarian assignment search was run by this verifier. No production or frozen diagnostic source was changed. Board tracking, this verifier directory, and requested run audit copies were updated; no commit was made.

## Run399 one-shot execution outcome

After run399 was reserved and linked, I invoked the exact command in its PLAN once. It exited with code 1 before any matrix evaluation because the manifest's source hash for `PREFIT_POINT_AUDIT.json` was stale: expected `9ac0678de91ca3eeed3ab799a1a1267ed96cf498302259e4694f442338c9e7f4`, actual `6a992a55d733e1b43be2c9adbaf8a2125e98f48a9b7ef8842f9022ccd677750a`. No `fd_result.json` was created; the script performed zero matrix calls, and no FD values or derivative results exist. In accordance with the one-shot reservation, I did not retry or alter its frozen inputs. Board101 was closed inconclusive and the execution record was copied to `experiments/runs/399/FD399_EXECUTION_AUDIT.json`. Any further FD check requires regenerated bindings and a new reservation/board.
