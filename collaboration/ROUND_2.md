# Second coordinated research round — 2026-09-26

## Reconciliation and ownership

Starting revision: `f1dcda6`; initial worktree clean. The handoff was read before
assignment. The current agent registry contained only `/root`; historical
handles were not inherited. Process command lines and `/proc` working-directory
inspection showed MCP services but no circuit optimizer or certificate job.
The ledger contained 371 runs, none running, and all 16 board hypotheses were
terminal. Initial ledger audit: 371 runs, 186 candidates, no problems. Initial
board audit: four submissions, hash/JSON checks only.

Completed exact13 work (runs370–371) was not scheduled again. Exact13 acceptance
review belongs exclusively to `verifier_round2`. Separate algebraic and numerical
assignments target changed twelve-CNOT constructions. Per-agent directories and
board events preserve ownership; the coordinator owns records and the commit.

The user interrupted then resumed this same round. The three new handles were
interrupted and explicitly resumed, rather than replaced. A process/ledger check
found no surviving optimizer or certificate job. The pending test session was
collected successfully. See `work/root_round2/RECONCILIATION.json` for provenance.

## Bounded assignments

- Algebraic: inspect canonical13 blocks for analytic cancellation or a distinct
  construction; structural investigation, no optimizer or repeated certification.
- Numerical: one twelve-CNOT topology changing at least two retained pairs from
  a single-deletion source; two seeded starts, at most200 iterations each,
  total compute cap120 seconds and one compute thread. Freeze and reserve first.
- Verifier: review completed exact13 evidence, branch/scalar/hash binding and
  native regrouping under explicitly separate gate models. No repeated execution
  of the existing full-matrix verifier.

## Coordinator checks

`python3 -m unittest test_collaboration_board.py test_check_circuit.py
 test_exact_check.py test_exact14.py test_experiment_log.py`: 15 tests passed.
`python3 -m unittest test_native_gate_check.py`: seven tests passed.
`python3 check_circuit.py collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`:
13 CNOTs, numerical maximum off-diagonal error `5.65198076432888e-16`, unitarity
error `4.440892098500626e-16`. This numerical check is diagnostic, not an exact
certificate. The gate list does not diagonalize total spin.

Coordinator checked the source hash against both existing exact result files:
`d500ff2a79c4717ff38d935d71b67bd798b06a7c6899a275fec25570b86b8d13`.
Both record256 empty exact residual dictionaries and60 local normalizations.
This inspection is not a new independent full-matrix implementation.

## Outcomes

All three researchers completed; board17–19 and ledger372 are terminal.
Reports are in `work/algebraic_round2/REPORT.md`,
`work/verifier_round2/REPORT.md` and `work/numerical_round2/REPORT.md`.

Algebraic work derived a fixed-angle terminal-CX deletion obstruction and the
bare chronological identity C01 C12 C01 = C12 C02. Intervening local rotations
prevent applying that identity directly to the current circuit. No smaller
construction was established. During adversarial review the coordinator fixed
the original column-support argument: Hermitian K gives nonzero K00 andK08,
so a hypothetical diagonal E with EK=KD would require E00=-1 and+1.
The researcher confirmed the corrected proof. This excludes only unchanged-angle
terminal deletion, with free output-label ordering, not refitted12CX circuits.

Run372 froze one12CX schedule: delete original slot6 and change original slots7
and8 from32,12 to02,03. Both seeded full-local fits fail independent gate-list
checks. Seed926201 reaches200 iterations with loss0.3352359084631581 and
maxoff0.7561025755556912. Seed926202 converges at110 iterations to
loss0.3125000000000019, maxoff0.999999999999997. Actual8.052 seconds stays within
120 seconds; no additional fits. One topology/two starts are bounded evidence,
not a topology exclusion or a lower bound.

## Accepted incumbent and gate-model costs

The independent verifier recommends accepting the preserved exact13 evidence;
coordinator acceptance is `EXACT13_ACCEPTANCE.json`. The exact prescription,
source hash and nonzero local/projective scalar transfer were reviewed. No
new exact certificate computation was launched. Runs370–371 reproduce the
same full-matrix verifier; separate relation/sign/scalar audits and100-digit
simulation are independent corroboration. The researcher's review used the
handoff directly and coordinator-transmitted source excerpts/audit summaries,
because its command launcher was blocked. No second exact matrix implementation
is claimed.

| Elementary two-qubit gate set (arbitrary one-qubit gates free) | Native count | CNOT compilation upper bound | Evidence |
|---|---|---|---|
|Fixed directedCX, all-to-all|13|13|Accepted exact radical prescription|
|CX or F(a,b)=exp(i a XX+i b YY), independently tunable a,b|10|14|Preserved exact14-derived baseline; native JSON numerical match and analytic provenance|

The hybrid baseline has sixCX plus fourF-family gates and dependency depth8.
It has no direct exact native-JSON certificate and is not a ten-CX circuit.
Its current numerical checker gives maxoff4.0149e-16 and phase-aligned compiled
matrix error8.4552e-16. No improved exact native regrouping was established.
The global CX interval is6<=Cmin<=13; six is inherited from earlier exact proof,
not re-proved this round. No optimality, total-spin, strong Schur or novelty claim.

README, SEARCH_14_STATUS, AGENTS and collaboration current records are updated;
old chronological search observations and exact14 certificates are preserved.

## Recommended next bounded round — not dispatched

Algebraic researcher: derive a locally decorated three-wire replacement using
C01 C12 C01 = C12 C02, tracking intervening rotations and available eigenspace
gauge. One structural investigation, with at most one reserved60-second exact
identity check if needed. Numerical researcher: only after a distinct analytic
plan, freeze two changed12CX schedules, one seeded start per schedule,
max300 iterations and120 seconds each (240 seconds total), one thread, one ledger
reservation. Verifier: independently inspect any passing saved gate list and
reserve a separate exact-certification budget before certification. No follow-up
is dispatched and no completed certification assignment is reused.

## Final checks and review

The22 focused tests passed. Coordinator reconstructed both saved numerical lists
independently, reran shared checker (expected exit1 for invalid candidates), and
checked source/code/snapshot hash binding. Native baseline checker passed.
Final ledger audit:372 runs,187 archived candidates, no problems; zero running.
Final board audit:four submission hashes/JSONs valid, no proof certification.
`git diff --check` passed. Adversarial review fixed the algebraic support argument;
no remaining confirmed defect. Remaining risks: exact13 uses one full-matrix
implementation reproduced, native JSON is numerical, bounded fitting failures
are not exclusions, and novelty scan is scoped to source deletion schedules.

All research JSON, scripts, reports, candidate snapshots and databases are kept
artifacts. Python caches/SQLite journals are ignored. No destructive cleanup,
remote action or push. Local commit follows final review; working-tree state is
reported at closeout.
