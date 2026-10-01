# Short-fit topology tournament from the fresh chain13 source

Ledger run 347 enumerated 65 one-slot topology proposals: each of 13 CNOT
slots replaced by each of the five other unordered wire pairs, with
ascending direction at even slots and reversed direction at odd slots.
It excluded exactly the 13 unordered schedules already fitted in run 340's
closing-edge experiment. This guard covers only
`search13_fresh_chain_closing_edge_result.json`, not every historical fit.
All 52 remaining proposals were frozen before fitting in
`search13_fresh_chain_fitted_tournament_frozen.json`.

Each short fit started from the run-339 source's local angles and varied
all 168 local parameters for at most 40 L-BFGS-B iterations. The six
distinct schedules with lowest **fitted** direct Frobenius loss were then
continued for at most 400 iterations. No ranking used raw mutation loss.
The source loss roundtrip error is `3.04e-18`, and there are no inherited
fixed local gates or eigenlabel restrictions.

| Proposal | Changed slot | Short loss | Deep loss | Deep maximum entry |
|---|---:|---:|---:|---:|
| 56 | 11 | 0.01896273041 | 0.01451941897 | 0.12483925699 |
| 39 | 7 | 0.07376410427 | 0.07322330470 | 0.35709414497 |
| 49 | 9 | 0.08996115480 | 0.08985588503 | 0.48004253362 |
| 47 | 9 | 0.08996785011 | 0.08985588503 | 0.48004253419 |
| 0 | 0 | 0.09703616361 | 0.00145526662 | 0.03519943734 |
| 30 | 6 | 0.12500086576 | 0.12500000000 | 0.70710678119 |

All 52 short candidates and six deep candidates independently failed
`check_circuit.py`, each with 13 CNOTs. Selection ranking and all frozen
schedules match saved gate lists. Saved-gate NumPy losses reproduce records
with maximum error `3.55e-15`. The full validation record is
`search13_fresh_chain_fitted_tournament_validation.json`.

Proposal 0 (`01->02` at slot 0) converged with a smaller maximum entry than
the source, but its Frobenius loss is worse than the source's
`0.00121583422757`. Proposal 56 reached its 400-iteration cap, also at a
worse loss. No ordinary-loss improvement or valid diagonalizer appeared,
so no further refinement was launched. A 40-iteration tournament can rank
local fits imperfectly; one warm start and one chosen direction per pair
do not exclude a topology or a circuit below 14 CNOTs.

Recorded command:

```bash
python3 search13_fresh_chain_fitted_tournament.py --short-maxiter 40 --deep-maxiter 400 --deep-count 6
python3 check_circuit.py search13_fresh_chain_fitted_tournament_p0_deep.json
python3 experiment_log.py show 347
python3 experiment_log.py audit
```

Preserve the saved result and frozen proposals when reproducing in a
separate workspace copy.
