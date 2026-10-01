# Architecture experiment with confirmed source defect

**Confirmed provenance defect:** intended source was run377 `curvature/base.json`.
Frozen config actually points to run376 `joint/target0.json`, SHA
`e8f7cf6193b6e466b6e62a22ab7f63331ebd30b82d72c6bf4dfe63df541dabae`.
The config was copied from curvature's INPUT configuration; its source was not
updated to curvature's OUTPUT. Ledger parent377 and selection prose therefore
misidentify the actual source. Independent verifier caught the discrepancy.
Frozen config, archived code, outputs and hashes are preserved unchanged.
Ledger finish notes and board explicitly record the error. The intended
run377-source assignment remains unexecuted; no rerun without a new budget.

Actual board40/run379: two chronological12CX unordered-pair mutations from
run376 joint0. Slot6 01→03 or slot8 12→02, I slot1 remains optimizer-only.
All168local coordinates free, saved56U3+12CX, q0 MSB/right cycle chronological,
ancilla-free all-to-all directed CX plus arbitrary locals. Independently reset
both at actual376 source; max200L-BFGS iterations each,120seconds total,
oneTorch/BLASthread, no random seed, labels, Hessian or SVD.

Historical read-only duplicate audit before reservation found no matches
among378ledger configs and explicit sequence arrays/gates in1394repository
JSON artifacts,13939sequence observations ignoring directions. Current
own/verifier architecture proposal files were excluded after a concurrent
proposal self-match was detected. Scope does not reconstruct implicit
generator rules or establish universal historical novelty.

Actual runtime6.051061241seconds. Slot6 mutation converged148iterations,
156evaluations, loss0.4999272375430227/maxoff0.5115281228445528. Slot8
converged49iterations/68evaluations, loss0.49999999999999956/maxoff
0.8904639939651671. Both12CX gate lists independently invalid. No reduction,
topology exclusion or result from the intended377 source follows.

Adversarial review confirmed the source-selection defect above. Source hash
binding accurately reflects actual376 input; it does not make the assignment
correct. No archive rewrite or hidden rerun. Other checks: full gate counts,
two genuine unordered-pair mutations, budget, axis checkpoints and independent
full-matrix agreement. Reports retain all evidence; Python cache ignored.
Coordinator owns commit. Correct next run requires explicitly selecting and
hashing the intended OUTPUT `curvature/base.json`, with fresh reservation and
new authorization/budget, rather than copying its input config source.

## Subsequent corrective outcome — 2026-09-26

After the stop state above, coordinator separately authorized corrective
run380. It used actual377 `curvature/base.json`, independently verified path/
hash before reservation/fitting, and fresh board44/ledger380. Both mutations
remain invalid; second optimizer ended `ABNORMAL: ` after48iterations without
restart. See `architecture_corrected/REPORT.md`. Thus the intended source377
assignment has now been executed under the new budget. Frozen379 code/config,
outputs, hashes and misleading historical parent metadata remain preserved.
