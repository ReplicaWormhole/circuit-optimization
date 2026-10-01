# Large-perturbation multistarts on two new two-slot topologies

Ledger run 336 selected run 335's cases 15 and 9 by their **final fitted**
Frobenius losses, respectively `0.06019279508950225` and
`0.10070215086680244`. Complete topologies were loaded from frozen
proposals. Each received two seeded starts from run 332's optimized source
local angles plus independent Gaussian noise with standard deviation 0.8.
The starts did not use run 335's fitted local angles. All 168 local SU(2)
parameters varied for at most 400 L-BFGS-B iterations per fit.

| Case | Seed | Final loss | Maximum off-diagonal entry | Iterations |
|---|---:|---:|---:|---:|
| 15 | 33700 | 0.45319808762535213 | 0.7170113996456816 | 210 |
| 15 | 33701 | 0.32841514786288330 | 0.7581633672793923 | 222 |
| 9 | 33800 | 0.40563166036033880 | 0.6199135756778068 | 400 |
| 9 | 33801 | 0.37500000000000266 | 0.9999999999999950 | 95 |

All four saved 13-CNOT candidates independently failed `check_circuit.py`.
Saved-circuit NumPy Frobenius losses match optimizer records within
`1e-12`. Case 9 with seed 33800 reached the iteration limit; the other
fits reported convergence. None improved the warm-start fits or found an
exact diagonalizer. Four starts are limited numerical evidence and do not
exclude the two topologies, other initializations, or a circuit below 14
CNOTs. The bounded experiment is complete.

Recorded command:

```bash
python3 search13_joint_prefix_two_slot_multistart.py --seed 33700 --maxiter 400
python3 check_circuit.py search13_joint_prefix_two_slot_multistart_case9_s0.json
python3 experiment_log.py show 336
python3 experiment_log.py audit
```

Preserve the existing result when reproducing in a separate workspace
copy. The script uses run 335's frozen proposals, so it does not regenerate
filesystem-dependent exclusion choices.
