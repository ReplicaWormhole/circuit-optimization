# Run 397: two full98 Hessians at a repeated low-loss point

Board 93, parent 396. Use exactly the best 98-coordinate vector of run 396,
seed 9261703. No draws, perturbations or fits. Target UV4=DU, qubit zero MSB,
chronological gates, arbitrary eigenvalue order and degenerate basis. Complete
word, parameter chart and dependency hashes are frozen in config.json.

At most two AD gradients and two 98-by-98 AD Hessians: original off-diagonal
objective and one fixed optimal spectrum-label branch. Compute the base raw
Hungarian assignment and sixteen constrained optima, each banning all slots
with one row's selected root. The minimum alternative cost minus base cost is
the distinct-root assignment gap. Same-root slot permutations do not count as
different branches. Fixed labels are held constant for derivatives. A strictly
positive gap permits interpreting that Hessian as local assignment-objective
curvature; rounding and perturbation limits remain.

One thread, 110-second internal / 120-second external cap, no optimizer.
Eigenvectors are stored as columns, with largest-magnitude entry positive.
Negative-curvature trigger: least eigenvalue below -1e-6. Any proposed direction
requires separately reserved independent finite-difference confirmation before
a separately reserved escape fit. No trigger is not a local-minimum certificate.
If the internal cap interrupts a Hessian, retain completed records and base list.

Code SHA: 2e518feeaa4a57bf6a3497a17696643049bf8123750f81628606478ac9f0ccda.
Config SHA: f70c638de7ac8fecd9ffc23750ae4b5b1ceecd92236855ecc3f31052f9933d1a.
Reservation and workspace precede independent preflight and derivative execution.
Save inputs, gate lists, losses, gradients, Hessians, eigenpairs, timing, outputs
and interpretation. Finish the ledger even if incomplete or unsuccessful.

```bash
timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/397/diagnostic.py > experiments/runs/397/stdout.txt 2> experiments/runs/397/stderr.txt
```

No numerical curvature observation provides an exact certificate, family
exclusion or incumbent change.
