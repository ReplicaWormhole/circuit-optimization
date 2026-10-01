# One repaired13CX anchor check

Board125. Preparation serialized one complete42gate primitive word, config and
hash manifest only; zero matrices/searches/checker invocations. No optimizer,
random draw, dependencies or shared checker changes. Seed0/no randomness.
The by-hand foursector review is ../REPAIRED13_STATIC.md.

Chronological CX allocation: Bellprefix2 + R2 + B=CSdag2 + P3 + exactCH1 +
finalCCH3 =13CX. P is S2,RCCX(3,0->2),Sdag2. Final CCH is
X1,Ry2(3pi/4),Rx2(pi/2),RCCX(1,3->2),Rx2(-pi/2),Ry2(-3pi/4),X1.
All rotations and phases use exact rational-pi strings. The elementary gate
list includes cx,h,x,s,sdg,ry,rx,u3 only, supported by both checkers.

Root must reserve the script hash with config_json containing both
`frozen_config` equal to config.json and `frozen_config_sha256` equal to its
byte SHA256 (additional reservation provenance keys are allowed). Copy
check_repaired.py,candidate.json,config.json,MANIFEST.json,PLAN.md unchanged
into experiments/runs/ID. Link actualID to board125 before dispatch.
The runner refuses a non-run directory, nonmatching/not-running ledger row,
script/checker/runtime/config/candidate mismatch or second invocation; the
exclusive execution_started.json marker is preserved on failure.

One numerical native_gate_check invocation has30seconds; one exact_check
invocation has60seconds. Each receives one-thread environment. The external
command cap is100seconds, allowing import/serialization overhead. Root runs
once only, preserving every stdout/stderr/result even on failure:

```
timeout --signal=TERM --kill-after=5s 100s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/ID/check_repaired.py > experiments/runs/ID/stdout.txt 2> experiments/runs/ID/stderr.txt
```

A pass must have numerical valid_diagonalizer=true and exact_diagonalizer=true
with successful checker exits, besides the frozen13CX count. It is a repaired
exact13 anchor for a prospective12CX compression/deletion experiment; it
changes no accepted incumbent and establishes no optimality or family exclusion.
Any12CX fit needs its own frozen reservation and resource budget.

Static adversarial review: literal chronology and all sector signs agree with
REPAIRED13_STATIC.md; config required/read keys are complete; source/candidate
hashes bind the actual complete word; timeouts and failed outputs are retained.
Full matrix verification remains unexecuted pending root reservation.
