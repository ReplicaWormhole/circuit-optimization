# First coordinated research round — 2026-09-26

## Infrastructure and repository

The project now has its own nested Git repository. Initial snapshot `6b7e033`
preserved the then-current project, including existing experiment work. Parent
history and tracking were left unchanged; no remote was configured or pushed.
One inherited trailing-space warning in a code snapshot was preserved.

`AGENTS.md` defines the team, evidence rules, compute limits and file ownership.
The coordinator implemented the transactional board; algebraic and native-gate
researchers worked concurrently and exchanged a useful decomposition directly.
The numerical researcher received a separate bounded construction hypothesis.

The board is a cooperative local CLI, not an authentication service or automatic
scheduler. Budgets are enforced by the coordinator and bounded search scripts.
Content hash audits check recorded artifacts, not mathematical correctness.

## Verified baseline

Coordinator reran `check_circuit.py`, `exact_check.py`,
`delete14_integer_audit.py`, and `algebraic14_raw_matchgate_certificate.py` on
the incumbent. Both independent gate-list exact checks establish `U V4 = D U`
with 14 CNOTs; the raw certificate establishes the four-block construction and
matches the gate list up to global phase. Numerical off-diagonal error was
`6.70e-16`. This is not an optimality or strong Schur-transform claim.

The native researcher saved an explicit hybrid CX+XX/YY baseline with six CX
and four independently tunable XX/YY gates: ten two-qubit gates, dependency
depth eight, constructive CNOT compilation upper bound fourteen. Coordinator
independently reran the native checker: off-diagonal error `4.01e-16`, matrix
difference from the compiled incumbent up to phase `8.46e-16`. Native JSON has
numerical validation and analytic provenance, not a direct native JSON exact
certificate. These counts belong to different elementary gate models.

## Algebraic handoff

`algebraic_hypotheses.md` proposes nine tunable XX/YY gates, nonmonomial changed
prefixes, and larger degenerate-sector gauges. It records the exact identity
`F(a,b)=E((a+b)/2) X_c E((a-b)/2) X_c`, explaining why unequal XX/YY blocks
cannot be silently counted as single equal-angle exchange gates. The native
researcher incorporated this distinction into its conventions and limitations.

## Checks and review

All 22 focused tests passed across collaboration ownership/dependencies/hash
integrity, native conventions, and existing circuit/exact/ledger functionality.
The ledger audit at its first coordinator snapshot checked 356 runs and 174
candidates with no problems. Scientific activity in another task continued
concurrently; this audit is timestamped evidence, not a guarantee about later
entries. Old run354 was subsequently finished by that task; we did not edit it.

Independent adversarial review of the implemented interface, policy, native
checker and records returned APPROVED with limitations. In reviewing numerical
preflight the coordinator found a wrong ledger-table name and required corrected
novelty provenance before reserving the run. Numerical outcomes are recorded below.

Unrelated concurrent research artifacts remain intentionally unstaged. The
shared scientific ledger and board are preserved as whole snapshots rather than
rewriting another task's records.

## Shared work arriving from another task

Another active task adopted this interface: hypothesis4/run359 produced
`collaboration/work/root/labelswap_svd_refine_candidate.json` with thirteen
CNOTs, and hypothesis5 assigned exact certification to its own collaborator.
Coordinator independently reran the numerical checker: it passes at maximum
off-diagonal error `2.00e-15`. This is numerical evidence, not an exact incumbent.
We did not duplicate that certification assignment. System `python3` could not
run the optional Qiskit crosscheck because Qiskit is not installed there; this
round makes no independently reproduced Qiskit claim for the new candidate.

## Numerical experiment360

Hypothesis3 froze two distinct complete 13-CNOT schedules, normal-random local
starts with seeds 9261303 and 9261304, full freedom in all fourteen local layers,
80 L-BFGS-B iterations per start, one compute thread, per-start cap55 seconds
within total execution cap120 seconds. Execution took3.77 seconds; both starts
reached80 iterations. The corrected novelty preflight read359 ledger configs
and1221 JSON files; the novelty scanner recognizes saved gate lists and singular `schedule`/`topology`
fields. It does not exhaust all configuration schemas, circuit equivalence
classes or unsaved trials. Superseded preflight is explicitly kept
under that name, not used as the experiment specification.

Neither candidate passes. Losses were0.1250000000604 and0.3790756296025; maximum
off-diagonal errors were0.7071067811631 and0.9950598775439. Coordinator separately
reconstructed the listed schedules and reran the checker; values agree with the
researcher's independent NumPy reconstruction. Run360 and hypothesis3 finished
inconclusive, with an explicitly invalid numerical diagnostic submitted and
handoffs recorded. This does not exclude either topology or establish any bound.

Final coordinator ledger audit checked361 runs and178 candidates with no
problems; board audit checked three submission hashes with no problems.
The planned first round is complete. Further searches need a new bounded
hypothesis; exact certification already underway in shared hypothesis5 belongs
to the other task. No unlimited follow-up loop or remote action was started.
