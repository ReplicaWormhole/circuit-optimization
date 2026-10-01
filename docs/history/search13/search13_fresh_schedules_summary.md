# Fresh full thirteen-CNOT schedules

Run337 tests two unordered interaction schedules without any inherited local chunks, six-CNOT prefix, or matchgate suffix. The first repeats the round-robin list 01,23,02,13,03,12 twice and appends 01. The second repeats 01,23,12 four times and appends 01. Each has two alternating CNOT-direction variants, making four directed schedules; arbitrary local gates make the direction variants equivalent in expressive power. There are two independent Gaussian local-axis starts per directed schedule (std0.8,2), seeds33800+10case+start, maxiter300. This amounts to four starts per unordered pair order, not four distinct architecture families. These starts are not Haar-distributed.

All eight saved candidates fail independent cycle verification. The chain case2/start1 improves the numerical frontier: direct offdiagonal Frobenius loss 0.00121583483759014, max offdiagonal entry 0.04303922739, stopping at the iteration limit. This is not a valid diagonalizer or a smaller verified upper bound.

Run339 refines that chain candidate from three starts (sigma0,.02,.1, seeds34000--34002, maxiter1200). All return to loss about0.00121583422757 and max offdiagonal about0.0430346; all are invalid. The source loss roundtrip is asserted, and all serialized candidate losses and CNOT counts were independently reproduced by NumPy and `check_circuit.evaluate`.

Review: all intended schedules are admissible under the necessary incidence/support-cone checks; each candidate has13CNOTs; all local layers vary; serialization preserves the fitted objective within1e-11; no exact or topology-exclusion claim. No confirmed defect found. Only the root-level retained search13*.json candidate files were used in the initial directed-topology duplicate guard. Replay in a fresh copy without these runs' own outputs; the complete hard-coded schedules and seeds are retained in the script/result. Both ledger code snapshots preserve the exact executed programs.

Commands:

```bash
python3 search13_fresh_schedules.py --seed 33800 --maxiter 300
python3 search13_fresh_chain_refine.py --seed 34000 --maxiter 1200
python3 experiment_log.py show 337
python3 experiment_log.py show 339
```

Scripts, candidates, result files, and this summary are intentional retained artifacts. The proposed next test adds 03 connectivity by rerouting one CNOT at each possible slot of the new chain source, retaining the13CNOT budget.
