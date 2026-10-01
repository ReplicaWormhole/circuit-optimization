# Independent curvature audit

Hypothesis34. No optimizer or certificate computation. The independent NumPy
axis matrix uses exp[-i(xX+yY+zZ)/2] with analytic sinc and direct bit-column
CX matrices; flat mapping is C-order(14,4,3). Chronological UV=DU, q0 MSB.
Source hash/chronology and source axis reconstruction upglobalphase1.19e-15
pass; polished saved-axis/fullgate agreement upglobalphase9.84e-16 passes.

Only base candidate was emitted:12CX, invalid, offdiagonal loss/16
0.009423003567534258, maxoff0.1523762940765612. Independent axis objective
0.009423003567534283 agrees. Least AD eigenvalue−3.343255264296238e-8 is above
trigger−1e-7; neither perturbation fit was run. Gradient norm2.0320596718e-7
is finite; no exact stationary-point assertion. Hessian symmetry error
3.33e-16, normalized vector error below1e-15, eigenpair residual4.85e-16.

Frozen independent central differences along saved least vector:
h=1e-3 gives curvature+4.748458570791314e-8; h=5e-4 gives
−1.3405943022348765e-8. All raw ±h objectives are retained. Differences are
around1e-14 and finite-step sign changes: this does not independently certify
negative curvature. Higher-order terms and roundoff can dominate this tiny
mode. Nearzero Hessian modes may be gauge or coordinate artifacts. PSD or
threshold-not-crossed behavior proves neither minimum nor topology exclusion.

Run377 observed1.865857seconds; guard120seconds between operations is not
preemptive during an individual Hessian. No extra search or FD step screening.
Run command: OPENBLAS_NUM_THREADS=1 python3
collaboration/work/verifier_label12/curvature/directional_check.py followed by
numerical config.json, curvature.json and verifier directional_check.json.
Fullgate check uses preceding verify.py and base.json. Checks pass numerically.

Adversarial review found no confirmed defect; tiny curvature sign is explicitly
unresolved. Exact13 remains accepted;12CX failure gives no upper bound.
Separate hybrid10native entangler baseline is not10CX. Intentional saved
scripts/JSON/reports/chats retained; coordinator owns commit, no push.
