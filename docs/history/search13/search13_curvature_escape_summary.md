# Curvature escape from run327

Run329 stopped before fitting because a finite-difference check of the largest Hessian eigenvalue used too large a step (1e-3). It was closed as inconclusive. Run331 repeats at step 1e-4 and checks both smallest and largest curvature modes against finite differences.

At the run327 source, the 168-variable label-free diagonalization loss is 0.0787330125161619 and gradient norm is 1.13e-7. The smallest computed Hessian eigenvalue is -3.67e-8; none is below -1e-6. The largest is 15.6424793. There are 117 eigenvalues below the positive-mode threshold 1e-4 (including small negative values); they are not established exact zero modes. This numerical calculation is consistent with a flat local minimum, but does not certify one or exclude a circuit at this topology.

Six fits started at ±pi times each of the first three Hessian eigenvectors with eigenvalues above 1e-4 (modes 117,118,119), with all local angles free and maxiter250. Every serialized candidate was rechecked independently with `check_circuit.evaluate`. All have 13 CNOTs and fail diagonalization. Best loss returns to 0.07873301251616045, max offdiagonal error 0.38665766296. Other starts end in distinct worse basins. No valid circuit was produced.

Reproduce in a fresh copy (the script guards against overwriting the result):

```bash
python3 search13_curvature_escape.py --maxiter 250
python3 check_circuit.py search13_curvature_escape_mode117_sign1.json
python3 experiment_log.py show 331
```

Review: source topology and loss roundtrip asserted; Hessian symmetry checked; extremal modes checked by central finite differences; fit losses checked by reevaluation; serialized candidates independently checked. Hessian and fits are floating-point evidence, not an exact certificate. No confirmed implementation defect remains. The rerun code and first-run code have separate immutable ledger snapshots.
