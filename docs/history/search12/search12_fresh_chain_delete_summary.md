# Twelve-CNOT deletion tests from the fresh chain lead

Run341 deletes the first or last CNOT from the run339 chain source. Each of the two12CNOT schedules is refitted from the source local layers and from an independent Gaussian perturbation of std0.3. Every local rotation varies, maxiter400; seeds34200,34201,34320,34321. The model retains14local layers, including two adjacent layers at the removed interaction; their product is another local layer, so this overparameterization adds no entangling resource.

All four candidates independently fail the original cycle checker and have exactly12CNOTs. Final direct losses range0.16265--0.41054. Best maximum offdiagonal entry is0.5688733941 (delete final interaction, perturbed start); it is not a valid circuit and does not improve the13CNOT numerical lead.

Checks: undeleted source loss roundtrip, exact CNOT deletion/count, optimizer loss reevaluation, independent serialized NumPy matrices and `check_circuit.evaluate`, agreement of saved losses <1e-11. Adversarial review found no confirmed defect. These four bounded fits exclude neither these topologies globally nor arbitrary12CNOT diagonalizers. The entire original target and eigenlabel freedom were retained.

Reproduce in a fresh copy without this result file:

```bash
python3 archive/legacy_root/search12_fresh_chain_delete.py --seed 34200 --maxiter 400
python3 experiment_log.py show 341
```

All scripts, candidates, results, and this summary are intentional retained artifacts.
