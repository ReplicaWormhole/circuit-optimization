# Deep refinement of the disjoint slot-8 topology lead

Ledger run 332 continued run 330's case 8, a 13-CNOT topology with slot 8
changed from `0->3` to `1->2` relative to run 327. The unchanged start and
seeded local-angle perturbations of standard deviation 0.02, 0.08, and 0.2
were fitted with all 168 local SU(2) parameters free, at most 1,000
L-BFGS-B iterations each. The source loss roundtrip was
`0.06019279537098438`, differing from the saved source optimizer loss by
`6.245004513516506e-17`.

| Perturbation | Iterations | Final Frobenius loss | Maximum off-diagonal entry |
|---|---:|---:|---:|
| 0 | 216 | 0.06019279505795965 | 0.25491354507909697 |
| 0.02 | 228 | 0.06019279505817157 | 0.25491354218619505 |
| 0.08 | 360 | 0.06019279505800743 | 0.25491356627468070 |
| 0.2 | 382 | 0.23501608509606697 | 0.88519561562618790 |

Every saved candidate independently failed `check_circuit.py`; all contain
13 CNOTs and have unitary errors below `2e-15`. The best maximum entry is
in `search13_joint_prefix_slot8_refine_trial1.json`. The first three trials
give evidence for a recurring invalid local basin on this topology. The
largest perturbation left it but reached a worse solution. None establishes
impossibility or an exact circuit below 14 CNOTs. Further improvement may
require a changed topology or a substantially different initialization.

Recorded command:

```bash
python3 search13_joint_prefix_slot8_refine.py --seed 33400 --maxiter 1000
python3 check_circuit.py search13_joint_prefix_slot8_refine_trial1.json
python3 experiment_log.py show 332
python3 experiment_log.py audit
```

The script refuses to overwrite the existing result file; preserve the
recorded output when reproducing in a separate copy of the workspace.
