# One fresh full98 fit

Reviewed on 2026-10-01. The user authorized the proposed fresh initialization
and one bounded fit. Initial worktree clean; ledger392 terminal, zero running.
GPT-6 Luna collaborators handled structural scope, numerical implementation
and independent verification; the coordinator reserved and executed the fit.
Hypotheses73–77 record the work. Proposal74 duplicated75 and was closed before
any computation. No prior research row, accepted circuit or checker was changed.

## Frozen search and independent preflight

Run393 freezes one vector from `np.random.default_rng(9261501).normal(0,.3,98)`.
All98 coordinates vary: 84 local rotation-vector coordinates, six retained
target-wire A1/B1 coordinates, and eight F angles. A1/B1 start nonidentity,
with norms0.43494048/0.13997853. This is a different coordinate point from the
run385/389 embedding; it does not guarantee a different gauge orbit or basin.

Chronology: `L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`.
F(a,b)=exp(i a XX+i b YY); qubit zero is most significant. Four CX plus four
F blocks compile to12CX. The objective is
`||offdiag(U V4 U†)||F² /16`, without fixed eigenvalue order or eigenspace basis.

One L-BFGS-B start used Torch automatic gradients, maxiter300/maxfun10000,
maxls30/ftol1e-15/gtol1e-10, one compute thread, an internal110-second deadline
and external120-second timeout. No restarts. The code hash is
`b19f1aaae55d705fb2487ecb04ca1f8e878688404c2e9a03133edb4e73b305aa`;
the full config hash is
`d347c2aeff25828c476354e10baeb0df184022eb8c1ecfc859e7bcc2e7b89b52`.
Reservation preceded independent preflight, which preceded execution.

The independent NumPy/SciPy initial coordinate/list matrix error was6.00e-16,
and native/compiled matrix error4.00e-16, up to global phase. The hash-bound
[prefit approval](work/verifier98fresh/PREFIT_APPROVAL.json) confirms the
literal seed, all98 free coordinates, gate counts and resource limits.

## Result: lower objective, invalid diagonalizer

| Quantity | Initial | Best recorded point |
| --- | --- | --- |
| Normalized Frobenius loss | 0.908501360688 | 0.454898570731 |
| Maximum target off-diagonal entry | 0.5295689784 | 0.6743630736 |
| Compiled CNOT count | 12 | 12 |

The fit returned normally after284 iterations and301 objective/gradient
evaluations. Measured optimization/serialization/check time was3.439s;
the external timeout bounded startup and the whole command. Best selection
includes all evaluated line-search points. The best gradient norm was3.90e-8.
SciPy's relative function-change criterion fired; this does not establish
stationarity, a local minimum or valid diagonalization. The lower Frobenius
objective coexists with a larger maximum off-diagonal entry.

Independent postfit checks cover all four complete initial/best native and
compiled gate lists, both vectors, target loss, unitarity, and source hashes.
Coordinate/list errors were at most1.08e-15; native/compiler errors at most
7.45e-16; unitarity errors at most1.56e-15. Independent best loss agrees within
5.56e-17. See [the audit](work/verifier98fresh/POSTFIT_AUDIT.json) and
[run393 conclusion](../experiments/runs/393/CONCLUSION.md).

No passing12-CNOT candidate was found. Run393 is terminal `failed` because the
bounded search found no diagonalizer; its optimizer execution succeeded.
The best full list is retained and archived in the ledger. No exact certificate
was attempted and the exact13 incumbent remains unchanged. Failure excludes
neither this family nor other12-CNOT constructions; global6<=Cmin<=13 remains.

## Closeout and remaining scope

No further fit is dispatched. A possible distinct next test is another
explicitly frozen nonzero-interior seed under a new finite budget. Repeating
unchanged seeds or launching an unlimited restart loop is not authorized.

The numerical script and config are unchanged after reservation. The independent
checker was extended with postfit checks; both its prefit and later code hashes
are recorded, and the later source is preserved as run393/verify_snapshot.py.
Scripts, inputs, four gate lists, history, receipts, reports, database updates,
candidate and code snapshot are intentional kept research artifacts.

Verification: the four focused unittest modules passed16 tests; ledger audit
passed393 runs/202 candidates; board audit found no integrity problems; static
compilation and diff checks passed. Adversarial review found draft logging and
hash-guard defects that were fixed before freezing/execution, and no remaining
confirmed defect affecting the run. Floating computation is not an exact
certificate or optimum proof. SQLite ResourceWarnings accompanied passing tests.
