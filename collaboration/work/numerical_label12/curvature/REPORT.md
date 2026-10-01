# Label-free refinement and curvature audit

Hypothesis35, ledger377. New source is run376 joint target0, frozen SHA in
config; all168local coordinates remain free on unchanged12CX slot1-deletion
schedule. Saved base has56U3+12CX, q0 MSB, chronological, ancilla-free directed
all-to-all CX and arbitrary locals. Shared checker and prior artifacts unchanged.

Reserved before all computation: at most100label-free L-BFGS iterations,
ONE168x168 automatic-differentiation Hessian of
`||offdiag(UV4U†)||F²/16`; only if lambda_min<-1e-7, two deterministic
±0.1Euclidean-normalized least-eigenvector starts, at most150iterations each.
Largest-magnitude eigenvector coordinate chosen positive. No randomness,
extra SVD, extra spectral targets or fallback. Torch/BLAS1thread,120seconds
overall cooperative cap including audit. Source/code/helper hashes frozen.

Runtime1.865857092seconds. Base converged47iterations/51evaluations with loss
0.009423003567534276 and maxoff0.1523762940765612: invalid. The Frobenius
objective improves slightly while the maximum entry increases; these are
distinct metrics. One AD Hessian produced:

- gradientnorm2.0320596718335954e-7;
- raw Hessian symmetry maxerror3.3306690738754696e-16;
- symmetrized least eigenvalue-3.343255264296238e-8;
- least eigenpair residual4.852236990128356e-16.

Threshold was not met, so ZERO escape fits ran. `curvature.json` retains full
point, gradient, raw Hessian, eigenvalues and normalized direction; result
retains full base axischeckpoint and gate-list hash. These are floating AD
values, not an exact mathematical Hessian certificate. Nonzero gradient and
tiny negative values do not establish a stationary local minimum or topology
exclusion. The assigned verifier independently checks full matrices and
directional finite differences with predeclared steps.

Adversarial review checked source/chronology/axis reshape, budget branch,
one-Hessian count, derivative-preserving Torch objective, eigenvector
normalization/sign, complete checkpoints and no unrequested targets. No
confirmed issue. Inherited rotation regularizer1e-24 avoids the zero-radius
AD singularity. Cooperative cap can complete one already-started Hessian
after its deadline; actual computation is far below120seconds. No further
search launched. Coordinator owns commit; all research output kept and
Python cache ignored. Recommend a distinct newly budgeted mechanism rather
than repeating this no-escape branch.
