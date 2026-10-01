# Frozen one-candidate analytic11CX check

Board137 proposed/claimed before preparation. Complete30gate rational-pi
candidate has11CX and onlycx,h,rz,ry,u3 primitives. Seed0,zero randomness,
no matrices/objective/checkers/optimizer run during preparation. The independent
fullphase derivation is SECTOR_REVIEW.md. Accepted12 remains preserved.

Root reserves code hash with config_json containing frozen_config equal to
config.json and frozen_config_sha256 equal to its bytehash; additional metadata
keys allowed. Copy check_analytic.py,candidate.json,config.json,MANIFEST.json,
PLAN.md,SECTOR_REVIEW.md unchanged into experiments/runs/ID; link actualID to
board137 ownerverifier11router. Runner checks activeledger status/code/config,
script/checker/runtime/candidate hashes, literal word and all deterministic
budgets before checker calls. Exclusive execution marker prohibits silent retry.

One native_gate_check invocation with60seconds,tolerance1e-9; one exact_check
invocation with90seconds; one-thread environment,external160second cap. Root
executes once only, preserving all outputs and failures:

```
timeout --signal=TERM --kill-after=5s 160s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/ID/check_analytic.py > experiments/runs/ID/stdout.txt 2> experiments/runs/ID/stderr.txt
```

Passing requires bothcheckerexit0 and validity flags, actual11CX and exact
16rowlabels[0,0,0,2,0,1,3,0,0,1,2,3,2,3,1,2] in roots(1,i,-1,-i).
Expected conductor32. The label guard is an independentprediction for this
specific word; arbitrary output ordering remains allowed in general. Acceptance
is root-owned and requires independent certificate/review; no optimality,
lowerbound, total-spin or familyexclusion claim follows. No unchanged accepted
source or sharedchecker was edited. This directory is intentional evidence.
