# Full coupled 98-coordinate diagnostic

Reviewed on 2026-10-01. The initial worktree was clean. Three collaborators
used GPT-6 Luna at the user's request: algebraic98 (structural map),
numerical98 (one bounded AD diagnostic), and verifier98 (independent matrices
and directional finite differences). Hypotheses 66–72 and scientific runs
389–392 retain the assignments, frozen plans, actual exchanges and outcomes.
No optimizer or incumbent certificate was rerun.

## Family and preflight

The exact control-wire absorption is a matrix-level equality of realizable
families, with free prefix boundaries remapped; it does not reduce CNOT count.
See [the derivation](work/algebraic98/REPORT.md). The tested chart has 84 local
rotation-vector coordinates, six interior target-q1 coordinates and eight F
angles. Its chronology is

`L0' CX01 A1 CX21 B1 CX31 L1' F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`.

Here F(a,b)=exp(i a XX+i b YY), qubit zero is most significant, and gate
lists are chronological. Four CX plus four F compile to 12 CX with arbitrary
one-qubit gates. The target is UV4=DU with free eigenvalue ordering and bases
inside degenerate eigenspaces.

The algebraic collaborator derived the boundary multiplication order. The
verifier independently checked a nonzero boundary-map instance and the
frozen zero-A/B embedding. Preflight caught and corrected serialization and
coordinate-conversion defects before reservation of the Hessian script.
The hash-bound [approval](work/verifier98/PREFIT_APPROVAL.json) checks the
actual corrected script: generic map error 4.61e-16, serialized matrix error
5.24e-16, and native-to-compiled error 3.86e-16, up to global phase. These
floating checks corroborate the structural derivation, not an exact circuit
certificate.

## Run 389: full coupled Hessian, no trigger

The ledger freezes seed9261101, the point inherited from run385 through
run387, the complete 98-coordinate chart, dependency hashes, method and
script SHA-256 `739c8585bcb5b0146a3d735dd06dc2c0be65a7eaf9c923ffffa766303e741957`.
Reservation and workspace creation preceded computation. The command used
one compute thread and an external 240-second timeout. The recorded 2.155s
timer covers loss/gradient/Hessian computation only; it stops before
eigendecomposition and serialization and is not total process runtime.

| Quantity | Result |
| --- | --- |
| Normalized off-diagonal Frobenius loss | 0.5922239996302492 |
| Full gradient norm | 1.60473627e-8 |
| Least Hessian eigenvalue | -2.69852057e-9 |
| Frozen negative-curvature trigger | lambda_min < -1e-6 |
| Maximum Hessian asymmetry | 3.33e-16 |
| Least eigenpair residual | 6.98e-16 |
| Eigenvalues with absolute value at most 1e-6 | 57 of 98 |

No eigenvalue met the trigger. Zero escape fits or restarts ran.
[Run 389](../experiments/runs/389/CONCLUSION.md) retains the input, complete
gradient/Hessian/eigensystem, script, config, stdout/stderr and full native and
compiled gate lists. Its sole diagnostic candidate is the inherited base
point, not an optimized improvement. Independent simulation confirms both
gate lists are unitary and equivalent; the compiled list has 12 CX but fails
diagonalization with maximum off-diagonal entry 0.9438619418. No valid new
candidate or exact certificate exists.

## Independent derivative check and preserved failures

Run390 stopped before any objective evaluation: the frozen input used the
first row of the eigenvector matrix, while its correct guard required the
first column. The original ledger note blamed the guard; the corrected
[conclusion](../experiments/runs/390/CONCLUSION.md) and board handoff identify
the input error without rewriting history. Run391 changed the guard to
accept that row and performed five evaluations on the wrong direction.
Its curvature near0.343 is not evidence about the least eigenvalue; that
attempt is explicitly failed and its [output](../experiments/runs/391/CONCLUSION.md)
is preserved.

Before the final reserved check, the coordinator independently confirmed
exact point equality, exact Q[:,0] equality, unit norm, source/code hashes and
the eigenpair residual. Run392 used precisely that frozen column with an
independent NumPy/SciPy objective, five evaluations, one thread and an external
120-second timeout. At h=1e-3 and5e-4, centered curvatures were respectively
2.06501483e-8 and8.88178420e-10. Both fail the -1e-6 trigger. These small values
do not resolve the sign of the near-flat curvature or verify all Hessian
entries. See [run392](../experiments/runs/392/CONCLUSION.md).

## Established scope and next hypothesis

This round completes the requested full coupled diagnostic at the frozen
point. It establishes neither exact stationarity nor a local minimum, family
exclusion, improved lower bound, or valid12-CNOT diagonalizer. The accepted
exact13 circuit remains unchanged; the established interval remains
6 <= C_min <= 13, with the optimum open.

The best next distinct hypothesis is a fresh, fully coupled98 initialization
with nonzero A1 and B1, rather than another same-point Hessian or unperturbed
run385 warm start. Proposed budget: one explicitly frozen seed, one
single-threaded fit, at most300 iterations and an external120-second cap,
no restarts. This is **not dispatched**. Freeze every initialization rule,
code hash and limit before reserving it; independently check any full gate
list and require an exact certificate before incumbent promotion.

## Review and retained evidence

Adversarial review followed `~/.codex/adversarial_review.md`. Confirmed
preflight defects were fixed before the Hessian run; the FD provenance defects
were preserved as failed attempts and corrected in a fresh reservation.
No confirmed remaining defect affects the frozen-point diagnostic. Residual
limitations include floating differentiation/near-flat signs and a remapping
helper chart corner at SU(2) -I that is not exercised at the frozen input.

The 23 repository tests and archive audit passed. Ledger and board audits
check provenance/integrity, not mathematical proofs. Scripts, JSON, reports,
run workspaces, candidate and code snapshots, and SQLite database updates are
intentional research evidence; Python caches are ignored. No dependency was
installed and no external state was changed.

Closeout independently compared the databases against the initial Git
snapshot: all388 historical attempt rows, all65 earlier hypotheses, and all359
earlier events are unchanged. Exact13 acceptance/source and exact14 baseline
are byte-identical. All four new workspaces contain plans and conclusions;
the ledger has zero running attempts and process inspection found no round14
scientific process. Latest child turn metadata confirms GPT-6 Luna for all
three collaborators.

Verification commands: the AGENTS.md seven-module `python3 -m unittest`
command (23 tests, pass), `python3 archive/verify.py` (1417 files/42 links,
no problems), `python3 experiment_log.py audit` (392 runs/201 candidates,
no problems), `python3 collaboration_board.py audit` (no integrity problems),
`python3 -m compileall -q` on the three owned work directories and four run
workspaces (pass), and `git diff --cached --check` (pass). The test suite emitted
SQLite ResourceWarnings while passing; these do not certify mathematical
correctness. No new scientific computation was performed during closeout.
