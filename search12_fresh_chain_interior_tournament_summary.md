# Interior deletions of the fresh chain13 lead

Ledger run 349 tests native 12-CNOT circuits formed by deleting source
slots 1 through 11. Endpoint deletions 0 and 12 were tested previously in
run 341 and are outside this run's scope. The complete source hash and all
native 12-CNOT schedules are saved in
`search12_fresh_chain_interior_tournament_plan.json`.

Optimization retains the source's 14 local layers and substitutes an
identity for the deleted CNOT, varying all 168 local coordinates. The saved
gate list omits that identity and contains exactly 12 CNOTs. Two adjacent
local layers at the deletion are redundant but do not use any additional
two-qubit resource. No fixed inherited local gates or eigenlabels are used.
The source loss roundtrip error is `3.04e-18`.

Each of the 11 deletions receives at most 80 L-BFGS-B iterations. The three
lowest **fitted** short losses, at deleted slots 1, 9, and 11, are then
continued for at most 500 iterations:

| Deleted slot | Short loss | Deep loss | Deep maximum off-diagonal entry |
|---|---:|---:|---:|
| 1 | 0.05024541756 | 0.05024047358 | 0.25881905382 |
| 9 | 0.08985599413 | 0.08985588503 | 0.48004255982 |
| 11 | 0.15118811881 | 0.15118811881 | 0.69430533711 |

All 14 saved short/deep candidates independently failed `check_circuit.py`.
Native schedule counts and the selection ranking were verified. Saved-gate
NumPy losses reproduce optimizer records with maximum error `1.33e-15`.
The full check record is
`search12_fresh_chain_interior_tournament_validation.json`.

The best native 12-CNOT candidate is
`search12_fresh_chain_interior_slot1_deep.json`, with maximum entry
`0.2588190538193609`; it remains invalid. No exact circuit or improved
ordinary source loss appeared, so the bounded task stopped. One warm start
per deletion and three local continuations do not exclude any topology,
all 12-CNOT circuits, or a circuit below 14 CNOTs.

Recorded command:

```bash
python3 search12_fresh_chain_interior_tournament.py --short-maxiter 80 --deep-maxiter 500 --deep-count 3
python3 check_circuit.py search12_fresh_chain_interior_slot1_deep.json
python3 experiment_log.py show 349
python3 experiment_log.py audit
```

Preserve the saved result and plan when reproducing in a separate workspace
copy.
