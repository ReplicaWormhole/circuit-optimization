# Frozen independent deletion-refit preflight

Board133; verifier11router. Preparation only:0matrix/AD/RNG calls. Root owns
reservation, execution and finish. Inputs preserve the frozen numerical11delete
family, fit source/config and certified406 candidate as byte-identical copies.

Exactly12chronological CX deletions. Every case checks the collapsed144-vector
base against the independent deleted-source gate product up to global phase;
base and perturbedinitial serializers; frozen family NumPy matrices; initial
Torch matrix; Torch objective/gradient against independent exponential matrices
and deterministic directional central differences. Independent primitives and
rightshift construction use MSB qubit0 and chronological multiplication.
Rotation vectors use scipy.linalg.expm(-i r.sigma/2), independently of the
family's sine/cosine implementation. No repository circuit checker is imported.

Exactly132 fullword builds=12*(8source/base/initial/serializer/matrix builds
+1TorchADloss build+2FD builds),12autogradgradient evaluations,24FDloss
values. These are one finitepreflight, not an optimizer or numerical search.
Matrix tolerance1e-9; objective tolerance1e-10; FD step1e-6; gradient tolerance
1e-6+1e-5*max(abs(ADdirectional),abs(FDdirectional)). Directions are
cos((j+1)*(deletionordinal+1)/37),normalized, j0..143, no direction randomdraws.
Preflightseed0; initializer reproduces onePCG64(10012101) generator with12
successive normal(0,.025,size144) draws, all1728values allocated before cases.
Everybase/noise/initial/direction vector and diagnostic is retained in outputs.

60seconds internal alarm and external cap, one-thread environment before module
imports, Torch threads/interops1. Exactlyone invocation; execution_started.json
exclusive marker prevents silent retry. No optimizer calls, fitting, restarts,
parameter updates, diagonalizer claim or postfit audit. Separate postfit checks
need a new reservation.

Root reserves the raw config.json object exactly, using config.method verbatim,
codehash of preflight.py; link actualID to board133 ownedverifier11router.
Copy preflight.py,config.json,MANIFEST.json,fit_family.py,fit_source.py,
fit_config.json,source_candidate.json,PLAN.md unchanged into experiments/runs/ID.
The guard binds all input hashes, actual runtime module versions/dependencies,
active raw-config ledger reservation and active board owner/run link before
any circuit matrices or random draws. Root then executes once:

```
timeout --signal=TERM --kill-after=5s 60s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/ID/preflight.py > experiments/runs/ID/stdout.txt 2> experiments/runs/ID/stderr.txt
```

Static adversarial review checked source collapse chronology/global phases,
serializer branch handling, all144 coordinate freedom, exact source/topology
coverage, initializer draw allocation, objective/gradient convention, strict
raw-config hash guards, timeout failure preservation and finite counts. No
blocking source issue found; runtime consistency remains unverified until
reserved execution. A preflight pass certifies consistency only, not an11CX
candidate, optimum or topology exclusion. Accepted12 remains preserved.
