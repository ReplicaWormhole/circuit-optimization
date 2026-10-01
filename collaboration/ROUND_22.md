# Round 22: ten-CNOT hypotheses from the certified eleven-CNOT construction

The incumbent is the exact eleven-CNOT word in run 411, independently certified
in run 412. The interval remains `6 <= C_min <= 11`. This round investigates
whether joint final reflection/selector synthesis or full-local refitting after
one CX deletion yields a ten-CNOT diagonalizer. Neither bounded failure nor a
restricted suffix obstruction establishes a general ten-CNOT exclusion.

All three collaborators use GPT-6.1 Sol. Algebraic10 owns structural board 139;
numerical10 owns frozen diagnostic board 140; the existing numerical_sol agent
acts as independent verifier10 on board 141. Coordinator board 143 records the
joint hypothesis. Preflight board 142 was accidentally claimed by coordinator22
before its new ID was checked; its event history is retained and the coordinator
closed its unexecuted reservation 413 after the board correctly rejected the
same-agent ownership mismatch. A fresh verifier10-owned preflight replaces it
at identical scientific scope; no computation occurred in 413.

Before numerical execution, the diagnostic is frozen in
`work/numerical10delete/PLAN.md` and `config.json`: source run 411; eleven
chronological single-CX deletions, preserving every other direction and order;
eleven full four-wire SU2 local layers, 132 free rotation-vector coordinates;
PCG64 seed 10012201, one perturbation of scale 0.025 per case; L-BFGS-B with
Torch analytic gradients, maximum 80 iterations and 4000 calls per case,
8 seconds per case, 96-second optimizer deadline, 105-second internal and
115-second external limits, one CPU thread. No restarts or deepening.

The fit script SHA256 is
`37528ea8ba197eeac28b874faacd0272a9d4c02d8b8e56f212f8b16c5dad4ec6`;
the family SHA256 is
`9a5f4929242590230925ba9f3a7cd428da6143c61d7e196eade4203f6fd37731`.
The configuration additionally freezes source, topology, plan, acceptance,
checker, dependency and runtime provenance. Root reserves and snapshots each
actual computation before invoking it. Independent preflight precedes fitting;
independent endpoint audit follows fitting. Any complete promising word requires
independent gate-list validation and an exact certificate before promotion.

Preparation and by-hand analysis are not fictitious scientific runs. Run IDs,
outputs, best candidate, limitations and closeout will be appended below.

Reservation 414 exited during its ledger guard before scientific execution: the
new guard mistakenly used target `cycle`; the ledger stores `V4`. Its failed
workspace and earlier preparations are retained. Both fit and preflight are
refrozen with `V4`; no matrices, RNG, AD or fitting occurred in 413 or 414.

## Scientific outcomes and limitations

Corrected preflight415 passed121 complete unitary builds,11AD gradients and
22finite-difference losses over all11cases in1.216s. Maximum phase-adjusted
matrix discrepancy1.21e-15; maximum directional gradient discrepancy8.69e-10.
The sole fit416 then completed all11cases in6.171s,834iterations and938calls.
Best deletion7 removes the last merged CX32 and retains tenCX; loss
0.013675822916672686,maxoff0.2044593930001434. No case passes1e-9 and this is
not a near-exact hit. The complete best word is `experiments/runs/416/best_candidate.json`;
canonical ledger copy is `candidates/dbd11b2e209b98c217cf595a9f8178d2279924ed013b45db5edb9f404f00edc0.json`.
All eleven literal deleted words, collapsed base words/local matrices, noise,
initial/best132-vectors and complete words, history and diagnostics are retained.

Independent audit417 passed77 complete unitary builds and all11endpoint,
noise allocation, topology, word-count, history, hash, best-selection and metric
checks in0.299s. Its114-file input manifest and separate primitive/SciPy-expm
implementation were frozen before its one60s one-thread invocation. No fit,
AD, recognition or exactness check occurs in this endpoint audit. The run
workspaces retain input/code/config/plan/dependency hashes, commands, outputs
and conclusions. Standard ledger ingestion and audit also recompute numerical
candidate metrics as bookkeeping validation; they are not optimizer restarts.

No exact10 word was proposed; no exact certificate or incumbent promotion is
justified. Exact11 and all earlier baselines remain unchanged. These failed
finite fits establish no topology exclusion, local-minimum proof or global
lower bound. Original fixed98-word membership remains unresolved.

## Structural outcome and next hypothesis

`work/algebraic10/SUFFIX_OBSTRUCTION.md` derives the joint target suffix blocks
I,Z,-iY,(Z-X)/sqrt2. With its prefix fixed and u/v sector labels preserved, a
two-interaction target-only replacement satisfies a corner product law and
cannot mix just the required fourth block. This is a narrow by-hand argument,
independently reviewed; changed earlier blocks, mixed sector labels and general
output permutations lie outside it.

The distinct unexecuted next hypothesis in `work/algebraic10/NEXT_HYPOTHESIS.md`
is the chronological schedule `02,13,01,23 | 02,32,10,02,12,32`: move the
controlled-H interaction before the last x-controlled merged reflection.
The bare CNOT triangle `02,10,12` equals `10,02`; incumbent rotation wrappers
prevent applying that cancellation directly. Joint wrapper and sector-phase
synthesis is the open question. Equivalence of the reordered and previous
full-local families is unestablished. Any next computation requires a new
finite budget, owned board entry and frozen ledger reservation.

## Review and closeout

Adversarial review found and corrected the new target guard defect before any
scientific execution. A draft endpoint-config key mismatch was corrected before
freeze. All unsuccessful/unexecuted reservations are closed and preserved.
The remaining structural proof is restricted as stated and received independent
by-hand review; no formal proof-assistant check or broader exclusion is claimed.
Twenty focused board/ledger/checker/native tests passed (existing SQLite
ResourceWarnings were non-failing). Ledger and board audits pass. Historical
ledger rows/events and accepted circuit artifacts are preserved. New scripts,
JSON, logs, databases and reports are kept research evidence; Python caches and
SQLite journals are ignored. Only team-assigned files are staged; concurrent
animation commits remain intact. No dependency installation or remote operation.

All five new ledger attempts are terminal; all round22 owned hypotheses are
closed. Historical unclaimed proposals86,88,94 remain preserved.

Final staged-artifact review also caught a post-audit edit to the frozen run416
conclusion. Its exact audited bytes were restored; all114 manifest-bound inputs
were rehashed and matched before commit. Runtime results and ledger rows were
unaffected. Later outcomes live in417 and this round report.
