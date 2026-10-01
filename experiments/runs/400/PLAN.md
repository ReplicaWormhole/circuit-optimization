# Run 400: refrozen independent finite differences

Board 102 replaces guarded zero-evaluation run 399. Same exact run 398 point,
fixed labels and two least-eigenvector COLUMN 0 directions, fully manifest-bound.
All nine dependency hashes were checked by verifier and coordinator before
reservation. The audit dependency is now an immutable copy. Its reconstructed
metadata version exactly matches the historical expected SHA; prior files and
run 399 were preserved. No numeric inputs changed.

Exactly ten independent NumPy/SciPy coordinate matrices: two centers and eight
signed samples at h=0.01 and 0.001. Each direction uses its matching original
off-diagonal or fixed-label objective. No reassignment, AD, full Hessian, fit,
random draw or direction re-estimation. Record raw first/second differences,
comparison to source gradient dot direction/eigenvalue, unitary error and calls.
Cancellation and finite-step effects limit interpretation at 1e-9 curvature.

One thread, 30-second hard external timeout and ten-matrix evaluation cap.
One invocation only; guards bind script, manifest, source result, point, labels,
direction columns and matching reservation ID before any matrix evaluation.
Script SHA: d0c1f42c51cd46fc00f3671bf43daddc902dc57e867967acd2410369ea968706.
Manifest SHA: 7b395a9e7b92e15b8a4a704d537b05fb78b896f895eb7a691ab17a951475aea1.
Runtime config SHA: e09ee0e7480f86a4e3f298e2a2641dea3d3bcf812409698f0edaa2aac00addbe.

```bash
timeout --signal=TERM --kill-after=5s 30s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/400/fd_check.py > experiments/runs/400/stdout.txt 2> experiments/runs/400/stderr.txt
```

Close the ledger after the single invocation, retaining failures as well as
outputs. No optimizer candidate or exact certificate results from this check.
The base12-CX candidate remains invalid. No curvature below the frozen -1e-6
trigger means no trigger-based escape fit. No-trigger is not family exclusion.
