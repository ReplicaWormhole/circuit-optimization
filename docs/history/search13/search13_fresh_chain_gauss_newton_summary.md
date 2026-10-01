# Fixed-eigenlabel Gauss--Newton from the fresh chain lead

Run342 infers nearest eigenvalue labels from the diagonal of U V U† for run339 trial1. Multiplicities are (6,3,4,3) for (1,i,-1,-i), matching V. The nearest-label margin exceeds1. The labels are held fixed while every local gate varies; the equation U V = D U still allows arbitrary basis choice inside each eigenspace. Other labelings are not excluded by this local test.

Two trust-region least-squares fits (source and Gaussian noise std0.03, seeds34500,34501) use an automatic512-by168 real Jacobian and at most35 residual evaluations. A seeded directional finite-difference check agrees to1.84e-10. Both fits terminate at squared normalized residual about0.001216502001576; their maximum original cycle offdiagonal entry is about0.04299887. Both13CNOT candidates independently fail the original cycle checker. The small change from0.04303458 reflects a slightly different objective, not a valid circuit.

The saved matrices independently reproduce the fixed-label residuals within1e-11. No confirmed review defect found. Solves and stationarity measures are numerical evidence only, not a certified local minimum or topology exclusion. The exact14CNOT incumbent remains unchanged.

Reproduce in a separate copy without the existing result file:

```bash
python3 search13_fresh_chain_gauss_newton.py --seed 34500 --max-nfev 35
python3 experiment_log.py show 342
```

Script, source labels, result file, two candidates and this summary are retained artifacts.
