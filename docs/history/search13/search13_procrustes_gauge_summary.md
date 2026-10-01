# Continuous output eigenspace gauge by block Procrustes fitting

Ledger run 324 tests a seven-CNOT tail after the first six CNOTs of
`topology14_exact_matchgate_rational.json`. The two tail schedules are cases
0 and 5 selected by `choose(..., 6)` from the 1,428 schedules that survived
the earlier **identity-gauge** necessary rank screen. Those ranks do not
certify compatibility after the gauge changes.

Let `T` be the exact eight-CNOT tail and `W` a seven-CNOT trial tail. Output
rows of the established diagonalizer split into shift-eigenvalue sectors of
dimensions 6, 4, 3, and 3. Any block unitary `G` acting within these sectors
preserves diagonalization. For fixed `W`, minimizing

`||W - G T||_F^2`

over the entire block-unitary gauge has a closed-form solution. In each
sector, take an SVD of the corresponding principal block of `W T†`, say
`L Σ R†`, and set the block of `G` to `L R†`. The implementation computes
this solution numerically, then fits all 96 local axis-angle parameters of
`W` to the fixed gauged target with L-BFGS-B. It repeats this alternation
three times and finally fits the label-free diagonalization objective.

Four fits used the deletion warm start and one seeded random start for each
of the two schedules, with 60 process iterations per alternation and 120
direct iterations. The best final trial is
`search13_procrustes_case5_deletion_warm_direct.json`. Independent
`check_circuit.py` reports 13 CNOTs, unitary error `5.12e-16`, and maximum
cycle off-diagonal error `0.5336726945493274`; it is invalid. All four
trials and intermediate losses are recorded in
`search13_procrustes_gauge_fit_result.json`. The result is a bounded
numerical failure for these schedules and starts, with a fixed prefix. It
does not exclude other output gauges, schedules, prefixes, or 13-CNOT
circuits.

Reproduce with:

```bash
python3 search13_procrustes_gauge_fit.py --seed 32600 --outer 3 --process-maxiter 60 --direct-maxiter 120
python3 check_circuit.py search13_procrustes_case5_deletion_warm_direct.json
python3 experiment_log.py show 324
python3 experiment_log.py audit
```
