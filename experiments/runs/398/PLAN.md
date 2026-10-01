# Run 398: repaired curvature diagnostic

Board 98 replaces unexecuted run 397. Exact point: run 396 seed 9261703's best
98-coordinate vector. Full native word and chart are in the frozen config;
compiled count is 12 CNOTs. No perturbations, random draws or fits.

At most two first-order AD gradients and 196 second-derivative rows construct
two 98-by-98 Hessians: original off-diagonal loss and the fixed optimal spectral
label branch. Seventeen raw-cost Hungarian assignments compute the distinct-root
gap, banning all slots of the selected root for each alternative row. This gap
ignores duplicate-root slot permutations. Eigenvectors are columns; their
largest-magnitude entry is positive. Negative trigger: eigenvalue below -1e-6.

One thread, 110-second cooperative internal deadline checked at assignment and
Hessian-row boundaries, 120-second hard external timeout. Internal timeout in
either loop is caught and partial records are serialized with base gate lists.
An external kill can still prevent serialization. No optimizer or hidden fit.
Derivative directions require a separately reserved independent finite-difference
check before a separately reserved escape fit. No trigger proves no local minimum
or exclusion. Complete outputs, versions, timing and limitations are retained.

Code SHA: 5e67087e57f9b7d3a8915ac41339db7877b22df958a1c60e67896456cac9b086.
Config SHA: a5e6ecbf165e548321c5b178edb0c4d2941cfa72dfbce2179fd89b95acaeae13.
Reservation precedes independent preflight, then computation. Preserve run 397.

```bash
timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/398/diagnostic.py > experiments/runs/398/stdout.txt 2> experiments/runs/398/stderr.txt
```

No numerical derivative or gate-list check supplies an exact certificate.
The 13-CNOT incumbent remains accepted until a complete 12-CNOT circuit is
independently validated and exactly certified.
