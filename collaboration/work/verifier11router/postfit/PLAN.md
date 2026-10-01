# Independent all-endpoint audit of run409

Board134. Root owns fresh ledger reservation and execution. Preparation read
static saved outputs and froze hashes only:zero matrices/derivatives/RNG calls.
Run409 has12completed cases; no endpoint passed1e-9. Near endpoints7and9 are
numerical observations only and do not change accepted12.

Frozen input_manifest.json binds93source artifacts under experiments/runs/409:
sourcecode/family/config/39gate exact12source/topology/whole result/allnoise/
globalbest/executionmarker, plus all12case vectors/deletedsource/base/initial/
best gate lists/case result/initial metrics. Every artifact read by the audit
is hashed. Runtime guard requires run409 closed with matching code/config.
Other saved checkpoints/local matrices are retained but not used by this audit.

Exactly84fullword builds=7percase: deletedsource gate product, independent
exponential base/initial/best vector products, and base/initial/best savedgate
products. Reuse independently reviewed primitive/SciPyexpm methods from
preflight408. No existing matrix checker/family/Torch is imported. Three
endpoint metrics percase are computed from these existing matrices: normalized
offdiagonal loss, maximum offdiagonal entry, unitarity and fourth-root error.

Audit all144coordinate vectors, source deletion ordinal/edge/order, full59gate
48U3+11CX words, globalphase equivalence to source/base, serializer equivalence,
loss/checker agreement and actual11CX counts. Replay onePCG64(10012101) with
12successive144normal draws sigma.025, allocating all1728before checkingcases;
verify exactnoise and initial=base+noise. No other random draws; auditseed0.
History/best selection/counts/80iteration/4000call caps, allcase results and
globalbest bytes/hash selection are independently inspected. Matrix tolerance
1e-9; loss/metric comparison tolerances1e-10; actual validity threshold1e-9.

One invocation only,60seconds internal alarm andexternalcap,one-thread env.
No optimizer, gradient, finite differences, parameter update, restart, exact
certificate or additional search. An auditpass verifies the saved numerical
evidence and provenance; it doesnotturn a nearendpoint into an exact11CX word.
Any followuprefinement/exactcheck needs a separate stated budget/reservation.

Root reserves the raw config.json object using config.method exactly and
script hash; link actualID to board134 ownedverifier11router. Copy audit.py,
config.json,input_manifest.json,MANIFEST.json,PLAN.md unchanged into
experiments/runs/ID. Source409 remains in place read-only. Root executes once:

```
timeout --signal=TERM --kill-after=5s 60s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/ID/audit.py > experiments/runs/ID/stdout.txt 2> experiments/runs/ID/stderr.txt
```

Exclusive execution marker forbids retries. All runtime/source/inputmanifest
hashes, closed409provenance and active raw-configledger/board link are guarded
before any noise or matrices. Partialresults/failure files retain completed
cases if interrupted. Static adversarial review checked all endpoint bindings,
noiseallocation/globalphase conventions, losses/validity distinctions, source
path safety, count/budget and failure handling. No blocking issue found;
runtime84matrixaudit remains unexecuted until root dispatch.
