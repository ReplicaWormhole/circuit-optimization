# Numerical collaborator round 1

Experiment 360; hypothesis 3. Two fresh complete 13-CNOT schedules with unrestricted one-qubit SU(2) layers. Fixed seed 9261303 (second start 9261304), 80 L-BFGS-B iterations per start, one thread, 120-second execution cap after imports and at most 55 seconds per fit. Frozen proposal is `plan.json`; `preflight_superseded_plan.json` explicitly records an earlier preflight that accidentally scanned the wrong ledger table. The corrected preflight scanned 359 ledger configurations and 1,221 saved JSON files, identifying 341 distinct ordered 13-CNOT schedules. Novelty is exact ordered schedule novelty in persisted records, not gate-equivalence novelty or absence from unsaved optimizer trials.

Both fits ended at their 80-iteration limit, without reaching a valid diagonalizer. Normalized offdiagonal squared norms were 0.12500000006042622 and 0.3790756296025249. Best maximum offdiagonal entry is 0.7071067811631252. These failures neither exclude the schedules nor establish a lower bound. No exact claim or incumbent change is made.

Commands: `python3 collaboration/work/numerical/fresh_nonmonomial.py prepare` (once before reservation); `python3 collaboration/work/numerical/fresh_nonmonomial.py run`. Reservation config and frozen source live in experiment 360's ledger/code snapshot. Re-running should use a new experiment/directory rather than overwrite these evidence files.

Verification: both serialized JSON gate lists were read back and evaluated through `check_circuit.evaluate`; a separate NumPy reconstruction directly rebuilds U3 gates, bitwise CX permutation matrices, and right shift and agrees with checker maximum errors and Torch optimization losses to tolerance 1e-10. Unitarity errors are 1.33e-15. This verifies serialization and numerical residuals, not mathematical correctness of an exact construction. No passing candidate exists here.

Adversarial review using `/home/lev/.codex/adversarial_review.md`: incorrect ledger table in preflight was fixed before reservation/run. No remaining confirmed issue in the bounded experiment; remaining limitations are two random starts, local-minimum dependence, imported optimizer code reused rather than independently rederived, and exact-order rather than equivalence-class novelty checking. The coordinator independently reviews the result and owns commits.

Cleanup: scripts, plans (including labeled superseded preflight), gate lists, result JSON, log, this note, ledger candidate copy and immutable code snapshot are intentional retained evidence. No old experiment, shared code, or user artifact was modified by this collaborator. SQLite ledger and board rows were updated through their public interfaces.

Actual fitting/checking elapsed seconds: 3.7729355610208586

Imported helper file SHA256 values recorded at closeout:

```json
{
  "search13_adaptive_topology.py": "14973daa02e3f0e97a8a8c8cc4d2c05f906bc1f2dedec4b2893f2cf39cceb471",
  "delete14_search.py": "b12320ef0491d0fc953d76aa92c39025f5ea9c40dddfa1458894d18eef5fd044",
  "check_circuit.py": "b5ee62d3b83d3d82bd7a47006ab7697ffcb57d892ff17decface1e7bcf1921b0"
}
```
