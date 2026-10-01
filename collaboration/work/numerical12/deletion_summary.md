# Twelve-CNOT deletions from passing numerical13 (run 363)

Hypothesis 7, owner `numerical12`, tests all 13 single-CNOT deletions of run359's numerically passing source, `collaboration/work/root/labelswap_svd_refine_candidate.json`. Source SHA256: `f020145fbc1681d1d67f01e046b2833e30e64714100090470e81afae14c0770d`. This uses a different solution from the invalid run339 source tested in earlier deletion experiments; no previous output paths were reused.

All 13 directed native schedules and optimizer schedules were frozen before fitting in `deletion_frozen.json`. The optimizer keeps all 14 local layers (168 free local coordinates), replaces the deleted entangler by identity, and varies ordinary offdiagonal loss `||offdiag(U V4 U†)||F²/16`. Serialized gate lists omit that identity and have exactly 12 CNOTs plus 56 U3 gates. Adjacent local layers are redundant but cost no additional counted CNOT. No eigenvalue labels were fixed.

Each deletion used the source local coordinates, max100 L-BFGS-B iterations. The three lowest **fitted** short losses, slots `[1,12,5]`, were deepened from their short endpoints at max500 iterations. No random starts or seeds. One CPU thread. Source reconstructed optimizer loss is `9.383234516547568e-30`; independent source NumPy loss differs by `3.002270330097649e-30`.

| Deleted slot | Short loss | Deep loss | Deep max residual | Deep iterations |
|---|---:|---:|---:|---:|
| 1 | 0.0502404736554 | 0.0502404735808 | 0.258819036999 | 44 |
| 12 | 0.129163367722 | 0.129163367721 | 0.706034789673 | 26 |
| 5 | 0.135723306351 | 0.135723304703 | 0.353553398125 | 64 |

All **16 saved candidates** independently failed `check_circuit.py`, each with 12 CNOTs. NumPy saved-gate matrices reproduce recorded losses within `1.34e-15`; maximum residuals, full schedules, source hash, code snapshot and best-three ranking independently match. All three deep fits reported convergence before their caps. Best candidate: `deletion_slot1_deep.json`, loss `0.050240473580837414`, maximum residual `0.2588190369994409`. No valid numerical12 or exact12 diagonalizer was found. The ledger and board hypothesis are closed inconclusive.

Commands, from repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/numerical12/deletion_tournament.py --short-maxiter 100 --deep-maxiter 500 --deep-count 3
OPENBLAS_NUM_THREADS=1 python3 collaboration/work/numerical12/verify_deletions.py
python3 experiment_log.py show 363
python3 experiment_log.py audit
python3 collaboration_board.py audit
```

Replay the fitting command in a workspace copy without these existing outputs. The script intentionally guards against overwriting frozen plans/results. Checkpoints are progress evidence, not a process-liveness claim or an automatic resume mechanism.

Adversarial review of the computation, source roundtrip, representation, serialized resource counts, selection and validation found no confirmed issue. Residual risk: one warm start and finite budgets test only this deletion neighborhood; numerical failure does not exclude a topology or prove a lower bound. All scripts, frozen schedules, checkpoint, result, validation, summary and candidates are kept research artifacts. No commits, status-file edits or follow-up experiments were made.
