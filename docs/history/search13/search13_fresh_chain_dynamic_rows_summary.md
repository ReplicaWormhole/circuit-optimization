# Dynamic positive-row weighting of the fresh chain13 lead

Ledger run 343 tests three starts on the fixed run-339 chain13 topology:
unchanged source angles and independent Gaussian perturbations of standard
deviation 0.05 and 0.15, seeds 34400--34402. Each trial has three stages.
Before each stage, compute the current squared off-diagonal residual
`r_i` in each output row and set `w_i = 1 + 9 r_i/max(r)`. Hold these
weights fixed during at most 120 L-BFGS-B iterations. After three updates,
polish the ordinary cycle loss for at most 250 iterations. All 168 local
parameters remain free; there are no inherited fixed gates or fixed
eigenvalue labels.

Every weight lies in `[1,10]`, so each fixed weighted loss has the same
exact zeros as the ordinary off-diagonal cycle loss. Weighted losses from
different stages use different weights and are not directly comparable.
All nine complete vectors and the residuals used to construct them are
recorded in `search13_fresh_chain_dynamic_rows_result.json`, together
with source hash, schedule, seeds and optimizer diagnostics.

The source loss roundtrip error is `3.04e-18`. The all-one weighted loss
differs from the ordinary objective by `4.34e-19`, and their gradients
agree exactly in the recorded check. All nine weighted-stage candidates
and three polished candidates independently fail `check_circuit.py`; each
has 13 CNOTs. Saved-gate NumPy weighted and ordinary losses reproduce
optimizer records within `1e-11`.

| Start perturbation | Polished ordinary loss | Maximum off-diagonal entry |
|---|---:|---:|
| 0 | 0.0012158342275662570 | 0.04303461098977682 |
| 0.05 | 0.0012158342275622577 | 0.04303463587894632 |
| 0.15 | 0.0012158342275702030 | 0.04303461200410960 |

All ordinary polishes returned to the source basin. The 0.15 start's
second weighted stage transiently reduced the maximum entry to
`0.03305317918397369`, but its ordinary loss increased to
`0.0013605491022734332`. It is not a valid diagonalizer or an improvement
of the ordinary Frobenius loss. This is a bounded numerical failure,
not a topology exclusion, certified minimum or lower bound. The exact
14-CNOT circuit remains the validated incumbent.

Recorded command:

```bash
python3 search13_fresh_chain_dynamic_rows.py --seed 34400 --stages 3 --weighted-maxiter 120 --direct-maxiter 250
python3 check_circuit.py search13_fresh_chain_dynamic_rows_trial0_direct.json
python3 experiment_log.py show 343
python3 experiment_log.py audit
```

Preserve the saved result when reproducing in a separate workspace copy.
