# Run396: fresh8 full98 eigenlabel starts

Board87 (replaces unexecuted proposal86); parent395. Seeds9261700..9261707.
For indexi, default_rng(seed) draws90normal coordinates with sigma[0.3,1,2,3][i%4]
then8uniform[-pi/2,pi/2] F angles. All98 free; one fit/start, atmost8starts.
Complete native word/coordinate chart and dependency hashes are in config.json.
FourCX plus fourXXYY compile12CX; target UV4=DU with arbitrary eigenorder/bases.

Reuse run395's independently validated raw row assignment function, frozen as a
dependency: spectrum6(+1),4(-1),3(+i),3(-i), raw-cost Hungarian with pinnedSciPy.
Per-evaluation selectedlabels held constant for gradient, recomputed next call.
Every assignment is recorded; duplicate-root permutations do not changeD.
Store initialization, bestvector, label+original losses, native/compiled lists,
hashes, checker, counts and caps. Nearhit<1e-18 stops batch for verification,
not certification. Piecewise smoothness may affect L-BFGS-B convergence.

One-thread sequentialL-BFGS-B max1000iterations/max20000evalperstart/maxls40,
ftol1e-16/gtol1e-12; total110sinternal/120sexternal. No extra restarts.
SourceSHAc0e5eeed37086b68f588d291c9fdfed69a71043e8fb52646c264bbc7b24c0c47;
configSHA775abe1ab572124817e907db842df1c938a5f4dbb0f2998c0c1aa34fb3a9194c.
Reservation/workspace precedes independentpreflight, which precedes optimizer.
Close the ledger with result and limitationseven for failedbatch.

```bash
timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/396/batch.py > experiments/runs/396/stdout.jsonl 2> experiments/runs/396/stderr.txt
```

No incumbentchange without independentgatevalidation and exactcertificate.
Numericalfailure is not exclusion, lowerbound or local-minimum proof.

Pre-execution review notes: the20,000 cap counts optimizer objective-gradient
calls. Each retainedendpoint also receives one scalar original-loss diagnostic
and native/compiled gate-list checks; atmost8 such original-loss diagnostics,
within the same overall timecap. These are not extra fits. If the deadline
expires in the tiny interval after the loop's prestart check but before the
first evaluation, the stopreason is no_evaluated_point; earlier checkpoints
remain valid. Deadline at the loop check is overall_deadline_before_start.
