# Run 361: bounded Euler gauge recognition

Target: ancilla-free four-qubit right shift, `UV4=DU`, q0 most significant, chronological13-CNOT source run359. The CNOT schedule and output labels are unchanged. The exact14 incumbent uses paired XX/YY blocks and rational phases including pi/3; this fresh chain13 schedule does not immediately inherit that block derivation.

Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/algebraic13_new/gauge_snap.py`.

One bounded experiment completed all12 stages, with at most400 evaluations per stage and one CPU thread. Nine Euler coordinates per stage were prescribed nearest pi/8 values, chosen by linearized nullspace sensitivity;108/168 coordinates were prescribed. All12 fits passed the numerical residual threshold; actual evaluations23–43 per stage. Squared normalized equation residual3.19e-28. The independent gate-list checker reports13 CNOTs, maximum off-diagonal error2.239e-14 and eigenvalue-root error3.55e-15 for `gauge_candidate.json`. `result.json` records every selected coordinate, residual and evaluation count.

The fully rounded `rounded_candidate.json` is **not** a diagonalizer: numerical maxoff0.418, and `python3 exact_check.py collaboration/work/algebraic13_new/rounded_candidate.json` returns `exact_diagonalizer:false`, conductor32. No exact13 certificate was obtained. Prescribed coordinates are represented numerically in the partially fixed candidate; its residual remains numerical evidence, not an exact claim.

A notable surviving coordinate is inputq1 phi=-pi/3 within floating precision; most remaining coordinates do not recognize as small rational multiples of pi. This suggests pi/24 rather than pi/8 as a possible new grid, but this was not searched. An alternative next step is analytic extraction of the remaining noncyclotomic algebraic parameters. Neither avenue has been established.

Adversarial review: no unrelated files changed, original source preserved, no new dependencies, no incumbent promotion. Relevant checks were gate-list numerical validation and exact rejection of the rounded candidate. Remaining risk: numerical gauge manifold following gives no theorem about the existence or absence of rational-angle or exact13 solutions. Run361 and hypothesis5 are closed inconclusive; all work artifacts retained intentionally. Coordinator owns commits and shared state.
