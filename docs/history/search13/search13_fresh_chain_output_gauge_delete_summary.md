# Output permutations followed by earlier-CNOT deletion

Run346 appends one of three output CNOTs (0->2,1->3,0->3) to the fresh chain source, then deletes CNOT slot0,4 or8 from the original13 interactions. All nine resulting native schedules have13CNOTs. Every local gate is free through14 arbitrary SU(2) correction layers; fixed chunks contain only the original source one-qubit gates and do not impose a local-gate restriction.

An output CNOT permutes the computational basis. Consequently it preserves exact diagonalization and permutes the offdiagonal entries of an approximate diagonalizer; this is the reason to try these output-label gauges. Appending before deletion gives a14CNOT source whose maximum residual exactly matches the original. Deletion can destroy diagonalization; it is not assumed to preserve it.

All nine bounded fits (maxiter250, source-derived warm locals, no random perturbation) independently fail the original cycle checker. The best maximum offdiagonal entry is0.34965966; all are worse than the chain source. No valid circuit or topology exclusion resulted.

Checks: source count, output-gauge residual invariance, split/merge counts and order, free-layer zero start, loss reevaluation, all13CNOT candidate counts, independent serialized NumPy losses within1e-11, original cycle validity checks. Review found no confirmed defect. Floating numerical failures remain local evidence only. Complete source hash, schedules, fitted angles and optimizer records are retained.

```bash
python3 search13_fresh_chain_output_gauge_delete.py --maxiter 250
python3 experiment_log.py show 346
```

Reproduce in a separate copy without the existing result file. Script, nine candidates, result and this summary are kept research artifacts.
