# Run 399: independent central finite differences

Board 101, parent 398. Exact frozen point: run 396 seed 9261703's best98 vector.
Each objective uses its own least-eigenvector COLUMN 0 from run 398. Manifest
contains full point, fixed labels, direction vectors and hashes. No eigenvector
row selection, direction re-estimation, random draw or optimizer.

Exactly ten independent NumPy/SciPy coordinate-matrix evaluations: two centers
plus both signs at h=0.01 and 0.001 for each of two directions. Match the original
off-diagonal direction to that objective and the fixed-label branch direction
to that branch only. Compare central first/second differences to saved AD
gradient dot direction / eigenvalue. No assignment at perturbed points, AD,
full Hessian or fitting. Finite steps and floating cancellation limit agreement.

One thread, thirty-second hard external timeout, ten-matrix evaluation cap.
Script SHA: 9d7926acf863d64c4f405db04f9a4e437390d2f5ba1842d080a10e92afe70816.
Manifest SHA: ae2acd12cbcf4d374c8556e8e406f635b660178854819d29aac3e29563db9c9f.
Runtime config SHA: f8a958aa44274ae2c02f1e52d2a70b46165fee27d3a18e3b0956abd34d88d8b0.
Script guards both hashes, source result/audit/point/directions and fixed labels.
Reservation, workspace and board link precede execution by assigned verifier.
Outputs go to fd_result.json with every sample, formula, error and call count.
No candidate optimization or incumbent change; base12-CX candidate remains invalid.

```bash
timeout --signal=TERM --kill-after=5s 30s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/399/fd_check.py > experiments/runs/399/stdout.txt 2> experiments/runs/399/stderr.txt
```

No reliable curvature below -1e-6 means no trigger-based escape fit. These two
direction checks prove neither local minimality nor family exclusion.
