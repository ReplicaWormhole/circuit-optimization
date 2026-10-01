# Run 395: one dynamic eigenlabel fit

Board84, parent394. Warmstart is exactly the best98 vector from run394 row
seed9261620; no RNG calls, perturbations or restarts. All98 coordinates free.
Native word L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6
compiles4CX4XXYY to12CX. Qubit0 MSB; chronological gates; target UV4=DU,
arbitrary eigenorder and orthonormal bases inside degenerate eigenspaces.

Objective min_D ||UV4-DU||F²/16 with repeated slots6(+1),4(-1),3(+i),3(-i).
Raw row costs, unperturbed SciPy Hungarian assignment, fixed slots/input order
and pinned SciPy version provide a reproducible tied-optimum rule. Labels are
held fixed for each gradient; all selected columns/losses are saved. Assignment
switches make the objective piecewise smooth; nonsmooth termination is possible.
Original offdiag objective and full gate-list checker remain acceptance checks.

One L-BFGS-B fit: max1000 iterations/max20000 objective-gradient evaluations,
maxls40/ftol1e-16/gtol1e-12, internal110 seconds/external120 seconds, one thread.
Nearhit label-loss<1e-18 interrupts fit for independent verification, not success.
CodeSHA27a8ff0747dbf79f555fe4c4dec429cde1cc5530644e9b9715f6f2289e4ffa50;
configSHA7e1a8e0528f347b3e05e76d4cc9c617f264c0be7d0d9d1273b2c48530f7ac8c2.
Dependencies and parent-resultSHA are frozen in config and reserved ledger row.
Independent finite preflight must pass before optimization. Record endpoints,
inputs, outputs, bestpoint, limitations and terminal status even for failure.

Command from repository root:

```bash
timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/395/fit.py > experiments/runs/395/stdout.txt 2> experiments/runs/395/stderr.txt
```

No exact certificate is implied by floating-point convergence. Keep13-CX
incumbent until independent full-list validation and exact certification pass.
