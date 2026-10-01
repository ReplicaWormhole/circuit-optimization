# Diverse one-slot mutations of the changed-prefix near-hit

Ledger run 330 performed 17 whole-circuit fits from run 327's refined
13-CNOT candidate. At every one of the 13 CNOT slots, it replaced the edge
by the complementary disjoint wire pair, using alternating orientations.
It also reversed the directions at slots 0, 4, 8, and 12. These proposals
are independent of raw numerical loss and are distinct from all six
mutations fitted in run 328. Each fit varied all 168 local SU(2) parameters
for at most 250 L-BFGS-B iterations.

All 17 saved candidates independently failed `check_circuit.py`. Case 8,
changing CNOT slot 8 from `0->3` to `1->2`, improved the direct Frobenius
loss to `0.06019279537098432` and maximum cycle off-diagonal entry to
`0.2549165175930817`. It is still invalid and stopped at the iteration
limit. Its gate list is `search13_joint_prefix_slot_diversity_case8.json`.
The result file retains every proposal and checked residual.

The four direction-reversal fits also failed; reversal at the final slot
returned to the old loss `0.07873301251616353`. The new slot-8 lead supports
testing edge changes across all slots instead of ranking proposals solely
by raw warm-start loss. None of these numerical results establishes
impossibility for a topology or an exact circuit with fewer than 14 CNOTs.

Reproduce:

```bash
python3 search13_joint_prefix_slot_diversity.py --maxiter 250
python3 check_circuit.py search13_joint_prefix_slot_diversity_case8.json
python3 experiment_log.py show 330
python3 experiment_log.py audit
```
