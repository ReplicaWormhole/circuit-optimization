# Spectral-label bounded round

Hypothesis28, ledger375. Gate model: arbitrary single-qubit unitaries and
directed all-to-all CX, ancilla-free; q0 most significant, chronological gates,
right cycle. Source is the saved run363 slot1-deletion lead, hash
`9a9dd0e0b430e935f99f4d998fb25d5e1b48160ad16a6c46a27e2574c59c2398`.
The deleted identity is retained only in the optimizer; all saved lists have
12CX and56U3. Adjacent local layers retain all168 axis coordinates.

Before screening, `config.json` and `search.py` were reserved and snapshotted
in ledger375 and linked to board28. Nearest multiplicity-constrained Hungarian
assignment uses roots(1,i,-1,-i) with counts(6,3,4,3). Four distinct-label pair
swaps rank by descending |B_rs|²+|B_sr|², then lexicographic pair order. Each
target starts independently from the same source: max100 L-BFGS iterations,
max2000 function evaluations, then at most3 truncated-SVD Newton corrections,
cutoff1e-6 and8 halving line-search trials. Overall computation cap120seconds,
including screening; Torch/BLAS1 thread. No randomness.

Actual runtime14.848544887seconds. Pairs(2,6),(3,7),(8,13),(9,12), all100
iterations and3 Newton steps, Jacobian retained rank60. Final fixed-target
losses from read-only full serialized-gate diagnostics:

| Pair | ||UV-DU||F²/16 | max offdiag(UVU†) |
|---|---:|---:|
|2,6|0.05111126056639771|0.2588190451027492|
|3,7|0.05111126056641396|0.25881904511477694|
|8,13|0.2883009523150334|0.2794011825428628|
|9,12|0.28830095231503283|0.2794011768358336|

All four gate lists are invalid diagonalizers. No reduction, exact claim,
topology exclusion, or lower bound follows. Fixed-target loss is not the
ordinary offdiagonal loss. `optimizer.loss` precedes Newton; `diagnostics.json`
contains final gate-list values. A rejected final Newton trial is not applied.
Source reconstruction matches the source full matrix up to global phase to
1.4488835837920476e-15. The shared checker supplies numerical validation only;
the assigned verifier provides the separate independent matrix check.

Reproduce only with a newly reserved run and fresh output directory. Existing
evidence is intentionally immutable. Read-only diagnostics command:
`python3 collaboration/work/numerical_label12/diagnostics.py`.

Adversarial review: checked source binding, exact chronological12CX schedule,
four-target scope, multiplicity preservation, independent warm starts, full
168 freedoms, budget, rejected-trial handling and final serialization.
No confirmed defect found. Wall-clock guard is cooperative at evaluation/
correction boundaries, so one in-progress Jacobian/SVD can overrun a deadline;
this run ended far below its120second cap. Hungarian ties depend on the current
SciPy implementation; pair-score ties explicitly use lexicographic ordering.
Recommendation: a new structural mechanism or bounded topology change rather
than unchanged repetition of these four fixed-label targets.
