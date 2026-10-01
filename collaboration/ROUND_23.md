# Round23: controlled-H earlier inside the merged ten-CNOT word

The accepted exact11 circuit and its run411/412 certificates remain the
incumbent. Round22's eleven single-CX deletion fits failed; independent417
endpoint validation passed. This new round tests a changed ordered topology,
not a restart of a closed fit. It also investigates whether bare CNOT triangle
absorption survives the rotations and sector phases of a complete diagonalizer.
Artifact types: math, computation, code and workflow. Qubit0 MSB; gates
chronological; target UV4=DU; arbitrary eigenvalue order and degenerate bases.

Coordinator board151 records the hypothesis before scientific computation.
The three GPT-6.1 Sol collaborators have distinct structural, numerical and
independent-verification assignments and exchange actionable handoffs directly
and through the board. Each uses a new per-agent directory; old evidence and
accepted circuits remain intact.

## Frozen scope before computation

New topology:
`02,13,01,23 | 02,32,10,02,12,32`.
The last two noncommuting edges in the first part of the old suffix are
reordered relative to run416 deletion7. Initial base132 rotation-vector
coordinates are transplanted unchanged from that case's best132 vector.
Every local layer is free, including those inside the decoder/router prefix;
only the CNOT edge order is fixed in the numerical family. The exact prefix
is held fixed only in the restricted algebraic argument. The case input SHA256 is
`853d53241573363a5baf9d86f833cdf9ec5c6304cc0b2b2fa168466015864ef7`.
The source is invalid, not a near-exact or certified ten-CNOT construction.
Full-local family equivalence or inequivalence is unestablished.

One bounded diagnostic, two starts only, all132 real local SU2 coordinates
free in11layers4wires3axes. NumPy Generator PCG64 seed10012301 allocates two
consecutive132-coordinate normal draws, sigma0.025 then0.25, at startup;
each initial vector is base plus its assigned noise. L-BFGS-B with Torch
analytic gradients, max250iterations/10000calls/maxls20 per case,
ftol1e-15/gtol1e-10,15seconds per case,35second optimizer deadline,
45second internal and60second external limits, one CPU thread. No restarts,
ranking, deepening or unfrozen follow-up. Base/noise/initial/best vectors,
complete54-gate10-CX words, histories, metrics and failures are retained.

Before execution, the numerical script, family, topology rule, source files,
code/config/plan/dependency hashes and runtime versions are frozen in the
ledger and copied to its workspace. TargetV4; parent416; board/ledger owners
must agree. Independent finite primitive/SciPy-expm preflight precedes fitting;
a separately frozen finite endpoint audit follows. Any promising complete
word requires independent exact certification before incumbent promotion.
Root reserves and executes; preparation performs no matrices, RNG, AD or fits.

The by-hand algebraic branch may produce a complete exact word or a restricted
obstruction. Optimizer failure would establish neither a topology exclusion
nor a lower bound or local-minimum proof. The interval remains6<=Cmin<=11
until independent exact evidence justifies a change. Original fixed98-word
membership remains a separate unresolved question.

Actual freezes, run IDs, outputs, limitations and closeout will be appended.

Preflight418 failed before any unitary, RNG or AD evaluation because its whole
JSON comparison treated intentionally changed descriptive metadata as gate
data. Failure0points/0matrixbuilds0.533s is retained. A fresh frozen semantic
`n`/`gates` comparator and owned board entry replace this failed packet.

## Actual freezes and outcomes

The sole main fit script SHA256 is
`6faca4de3809b39bdaf69e29f544ae25f157af4e09e8e37008a53e93f8205fce`;
config SHA256 `055292df8e3481f97abd697dc13190ba838b8afd0557e83c7b8fb5073d18d083`.
The family SHA256 is
`9a5f4929242590230925ba9f3a7cd428da6143c61d7e196eade4203f6fd37731`.
Seed, topology, method, source, runtime/dependency, plan and configuration
provenance was reserved before execution. An earlier unreserved packet's
no-ranking wording was clarified to prohibit adaptive allocation while
allowing a deterministic global-best summary after both fixed-budget cases.
Superseded preparations are retained; they are not scientific runs.

Corrected preflight419 passed 22 complete unitary builds, 3 AD gradients and
6 finite-difference losses at the base and two initial points in 1.110 seconds.
Maximum phase discrepancy is 9.55e-16; directional gradient discrepancy 1.05e-9.
It independently reconstructs the literal edge swap, preserving all primitive
locals, and verifies raw source/base/noise allocation and serializers.

Fit420 completed both starts in 0.942 seconds: 136 iterations and 152 calls.
The unchanged raw vector has new-topology base loss 0.8144908421443564,
rather than the old topology's 0.0136758. Initial losses are 0.8040098930409487
and 0.8919066054116255. Case0 (sigma0.025) ends at loss 0.12500000000000047,
maxoff 0.7071067811865476; case1 (sigma0.25) at loss 0.2627760545065096,
maxoff 0.6361033566248393. Both fail 1e-9. Global best is selected by loss,
not maximum off-diagonal error: `experiments/runs/420/best_candidate.json`;
canonical copy `candidates/ff326b87764f17d85df55d737f37b7d39717682f060703a6a3df31a895a300a8.json`.
Every base/noise/initial/best vector and full 54-gate, 10-CX word, history,
metric, source and command is retained. No restart, adaptive allocation,
deepening or recognition was dispatched.

Independent audit421 passed 14 full unitary builds over both cases in
0.173786 seconds, one thread. Its 43-file manifest, 264-noise replay and separate
primitive/SciPy-expm implementation were frozen before execution. Full words,
raw132 vectors, semantic literal transplant, topology, histories, counts,
metrics and global-best selection agree. It explicitly distinguishes the
minimum-loss and minimum-maxoff cases; both are invalid. The immutable ledger
finish notes' 0.095-second elapsed figure is a clerical error; saved result/stdout
are authoritative and `experiments/runs/421/ERRATUM.md` corrects it. No numerical
or exact conclusion depends on that figure. No new exact candidate exists.

## Algebraic result and next hypothesis

The independently reviewed by-hand report
`work/algebraic10reorder/PHASE_OBSTRUCTION.md` retains the exact uv10 operator
`W10=[[0,I],[Z,0]]`: the lower Z contains the y-dependent quarter-turn phase.
With computational sector/control bases fixed, one controlled reflection has
relative target block E with scalar E^2. An early x-only CH could diagonalize
x only if E P and E-dagger were proportional, forcing the non-scalar
P=Q0 Z Q0-dagger to be scalar. Later x-block-diagonal gates cannot remove the
remaining x-off-diagonal block. This excludes only the stated template, not
the full132 family tested here. Coordinator and verifier checked phases,
commutators and invertibility by hand; no formal proof-assistant check ran.

The best next hypothesis is to solve W00=CZ and W10 jointly using tilted
x control axes and a later non-diagonal x readout rotation after the second
x-to-y interaction, then extend to both v=1 sectors using the two v-controlled
reflections. A phase-compatible v=0 construction would provide an informed
initialization. No solving angles, complete ten-CNOT word or new run budget is
asserted. Original fixed98 membership and general ten-CNOT existence remain
open. Finite failed fitting proves no local minimum, topology exclusion,
family equivalence/inequivalence or global lower bound. Exact11 remains accepted.

## Closeout

All four new attempts are terminal. The metadata comparator was corrected and
independently reviewed; failed418 remains preserved. Unused proposed board152
has an inconclusive dependency and cannot be claimed; it is administrative
history, not a pending computation. Twenty focused board/ledger/checker/native
tests passed, with existing non-failing SQLite ResourceWarnings. Independent
preflight/endpoint audits and the historical-row/hash checks provide numerical
and provenance evidence, distinct from an exact certificate or a broad theorem.
Research JSON, scripts, snapshots, logs and databases are retained intentionally;
Python caches/journals ignored. No dependencies or remote state were changed.

Final adversarial review found no unresolved blocker after correcting the
metadata comparator and documenting the elapsed-time erratum. A stray plus
in the new status prose was removed during diff review. All 417 previous
ledger rows, 147 previous hypothesis rows, 760 events and four submissions
remain byte-for-byte equal at the row level to the starting Git snapshot.
Accepted11/12/13/14 gate/certificate files and all 43 audit-bound inputs match
their original bytes. New research files are intentional retained evidence;
only this team's files are staged. Timeout-recovery branches were not exercised
by the successful short runs; the mathematical obstruction is restricted and
reviewed by hand, with broader existence and exact optimality still unresolved.
