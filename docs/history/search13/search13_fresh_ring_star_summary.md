# Fresh ring and star thirteen-CNOT fits

Run345 extends the fresh-schedule driver with an optional ring/star family. Default round-robin/chain behavior is retained; run337's exact earlier code is preserved in its immutable ledger snapshot. The ring repeats01,12,23,03 three times and appends01; the star repeats01,02,03 four times and appends01. Each has two alternating direction variants and two independently seeded local-axis Gaussian starts (std0.8,2), giving eight fits with maxiter400. These are two unordered interaction orders, not four distinct expressive architecture families. Arbitrary local gates can absorb changes of CNOT direction.

No inherited prefix, suffix or local chunks are used. All eight serialized candidates have13CNOTs and independently fail the original cycle checker. Losses range0.44093--0.51123; best maximum offdiagonal entry is0.4492480. One star fit terminates with an abnormal optimizer line-search status, recorded verbatim; its saved loss was independently reproduced, so this is a bounded failed fit, not a convergence claim or topology exclusion.

Checks: necessary incidence/support admissibility, limited duplicate guard against retained root-level search13 candidate JSONs, output schedules and counts, optimizer loss reevaluation, independent serialized NumPy losses within1e-11 and validity reevaluation. Adversarial review found no confirmed issue. Full schedules and seeds are in the result; all outputs are intentional retained research artifacts.

```bash
python3 search13_fresh_schedules.py --family ring-star --seed 34700 --maxiter 400
python3 experiment_log.py show 345
```

Reproduce in a separate copy without these runs' own output files; the driver refuses to overwrite results and excludes previously retained matching candidate schedules.
