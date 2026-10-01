# Opening02 tournament recovery (run 354, parent 351)

Run 351's original handle 1646 was unavailable and OS inspection found no process. Its frozen plan and seven saved short candidates (proposals 5–11) were preserved and independently checked. Run 351 was closed inconclusive due to interruption. Run 354 reused those endpoints without rewriting them. Their optimizer iterations, convergence flags and initial losses are unavailable; reconstructed losses are explicitly labeled as saved-gate evaluations.

Source: `search13_alternative_opening_refine_trial0.json`, ordinary loss `0.0014552666221479514`, roundtrip error `3.04e-18`. Its SHA256 is recorded in the result. All 168 local parameters were free, with 14 pure local layers and identity fixed chunks. New short fits used source angles; deep fits used their corresponding short endpoint. No random seeds or eigensector restrictions.

The frozen scope contains 65 one-slot other-pair proposals, five omitted and 60 retained. Directions ascend at even slots and reverse at odd slots. Exclusions cover only 66 distinct unordered schedules from run347's frozen proposals, omissions and original source, plus run340's saved closing-edge schedules. Full schedules, input hashes and exclusion records remain in `search13_opening02_fitted_tournament_frozen.json`. Every retained schedule changes two edges relative to the original run339 chain.

Run 354 performed **53 new short fits** at max40 iterations, reused **seven original short endpoints**, then deepened the **six lowest fitted direct losses** at max400 iterations. Ranking and frozen schedules independently match the saved gates.

| Proposal | Slot | Short loss | Deep loss | Deep maximum entry |
|---|---:|---:|---:|---:|
| 7 | 1 | 0.01483931469 | 0.01360587720 | 0.13439572601 |
| 56 | 11 | 0.02628861510 | 0.00911730094 | 0.10922791848 |
| 39 | 7 | 0.06607670729 | 0.06114697494 | 0.32021474825 |
| 43 | 8 | 0.07221398462 | 0.07221387488 | 0.37463658635 |
| 47 | 9 | 0.07902352343 | 0.07898964451 | 0.38951821268 |
| 34 | 6 | 0.10118848622 | 0.09461038033 | 0.28511643254 |

All **66 endpoints** (60 short and six deep, including **59 newly saved candidates**) independently failed `check_circuit.py`, each with 13 CNOTs. Saved-gate NumPy losses reproduce records within `2.67e-15`. Recovered hashes, code snapshot, ranking and every schedule were verified. Every deep fit reported convergence before its cap. The best candidate, `search13_opening02_recovery_p56_deep.json`, has loss `0.009117300938607778` and maximum entry `0.10922791847857725`, worse than source. No valid circuit or improved source basin appeared; the bounded task stopped.

## Commands and limits

```bash
python3 search13_opening02_recovery.py --short-maxiter 40 --deep-maxiter 400 --deep-count 6
python3 search13_opening02_recovery_verify.py
python3 experiment_log.py show 354
python3 experiment_log.py audit
```

For replay, use a workspace copy with the frozen inputs and seven original gate lists. The script writes an atomic checkpoint after each endpoint. This first recovery started without a recovery checkpoint and finished uninterrupted. **Future resume limitation:** checkpoint reuse does not bind newly saved endpoint hashes or CLI iteration limits to the original config. Validate those before any later resume. The current final verifier checks every schedule/loss, recovered original hashes, code snapshot and selection ranking.

Adversarial review found no defect affecting this completed first recovery; the checkpoint binding limit remains documented. One warm start, one direction per pair and finite iteration budgets do not exclude a topology or prove a CNOT lower bound. Candidates, checkpoint, script, verification script, result, validation and summary are intentionally kept research artifacts. No status file edits or commits were made by this agent.
