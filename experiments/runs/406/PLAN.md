# One explicit12CX candidate check

Board127 proposed/claimed before this candidate12 pack was created. A preliminary
merged12 serialization exists under the previous by-hand review scope; preserve
it and the unexecuted repaired13 pack. This pack is the new authoritative frozen
candidate12 input. No matrix, objective, optimizer, numerical checker or exact
checker has been run during preparation. Seed0, zero randomness.

Complete39gate primitive word has12CX: Bellprefix2 + router2 + mergedBP4 +
exactCH1 + repairedfinalCCH3. The BP merge uses controls(x=q0,v=q3), targety=q2:
Rz_y(-3pi/4),CX02,Rz_y(3pi/4),Ry_y(3pi/4),CX32,Ry_y(-3pi/4),
Rz_y(-pi),CX02,Rz_y(pi),Ry_y(3pi/4),CX32,Ry_y(-3pi/4),phase_x(pi/4).
The final CCH remains X1,Ry2(3pi/4),Rx2(pi/2),RCCX(1,3->2),
Rx2(-pi/2),Ry2(-3pi/4),X1. See ../MERGE4_STATIC.md and
../REPAIRED13_STATIC.md for independent by-hand identities/sector checks.

Root must reserve the script hash with config_json containing `frozen_config`
(equal to this config.json object) and `frozen_config_sha256` (its byte hash).
Additional reservation metadata keys are permitted. Copy check_merged.py,
candidate.json,config.json,MANIFEST.json,PLAN.md unchanged into
experiments/runs/ID and link the actualID to board127 before dispatch.
The script checks active ledger status, source/candidate/config/runtime hashes,
complete literal word equality and deterministic budgets before checker calls.
An exclusive execution_started.json marker prevents repeat execution, including
failed first invocations. Preserve all logs/results/failure evidence.

Exactly one native_gate_check invocation:30seconds,tolerance1e-9.
Exactly one exact_check invocation:60seconds. Single-thread environment for
both; external100seconds. Root executes once only:

```
timeout --signal=TERM --kill-after=5s 100s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/ID/check_merged.py > experiments/runs/ID/stdout.txt 2> experiments/runs/ID/stderr.txt
```

Both checker exits and validity flags must pass. The exact checker output labels
must additionally match the independent predicted row labels
[0,0,0,2,0,1,3,0,2,1,0,3,2,3,1,2], in roots(1,i,-1,-i), as frozen in config.
Arbitrary output eigenvalue order is allowed by the task; this label guard is
an independent consistency check on this particular construction, not an added
general diagonalizer requirement. Mathematical acceptance remains coordinator
owned and requires review of the complete saved word and exact certificate.
No optimality/lower-bound/family-exclusion claim follows from a passing12CX word.

Static adversarial review: checked all four block products, target/control
rotation chronology, relative T_x phase, all sector actions, primitive support,
count, required config/read keys and strict ledger/config binding. No unresolved
issue found before execution. Runtime/matrix verification is pending reservation.
