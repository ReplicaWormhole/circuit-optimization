# Positive output-row residual weights

Ledger run 338 tried four fixed row-weight vectors on each of two 13-CNOT
topologies: the run-332 source and frozen run-335 case 15. The vectors are
residual-adaptive (`1 + 9 r_i/max(r)`), reversed emphasis (`11-w_i`), and
two seeded exp-normal vectors clipped to `[1,10]`. All complete vectors,
seeds, source hashes, source schedules and row residuals are saved in
`search13_positive_row_weights_result.json`.

The weighted loss is `sum_i w_i sum_(j != i) |(U V4 U†)_ij|^2 / sum_i w_i`.
Every weight is strictly positive, so its zero set is exactly the ordinary
cycle off-diagonal loss's zero set. No eigenvalue labels are imposed.
All 168 local parameters vary in each stage: at most 200 weighted
L-BFGS-B iterations followed by 250 ordinary-loss iterations.

Both sources roundtrip correctly (errors `1.11e-16` and `4.16e-17`). The
all-one weighted objective and gradient also match the ordinary objective
within `1e-12`. All 16 saved stage candidates independently fail
`check_circuit.py`; each has 13 CNOTs. Saved-gate NumPy weighted and ordinary
losses match optimizer records within `1e-11`.

The adaptive weighted stage lowers the largest off-diagonal entry to
`0.24838542635924582`, but remains invalid. All ordinary polishes return to
the same loss near `0.06019279506` and maximum entry near `0.254914`. The
eight fixed weight vectors did not escape that ordinary-loss basin. This
is bounded numerical evidence, not a certified minimum or topology
exclusion. No exact diagonalizer below 14 CNOTs was found.

Recorded command:

```bash
python3 search13_positive_row_weights.py --seed 33900 --weighted-maxiter 200 --direct-maxiter 250
python3 check_circuit.py search13_positive_row_weights_t0_w0_direct.json
python3 experiment_log.py show 338
python3 experiment_log.py audit
```

Preserve the saved result when reproducing in a separate workspace copy.
