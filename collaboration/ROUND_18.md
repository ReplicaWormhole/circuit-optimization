# Round 18: plateau curvature and the next ordered family

Reviewed 2026-10-01 through run 400. The valid 12-CNOT target remains unmet;
the accepted exact 13-CNOT incumbent and exact14 baseline are unchanged.

## Completed evidence

Run397 was reserved and frozen but never executed: static review found that an
assignment-loop timeout could bypass output serialization. It is closed failed
with zero scientific evaluations. Replacement398 preserved those artifacts,
repaired the catch, and completed both full98 AD Hessians in 4.839 seconds,
within 110-second internal / 120-second external one-thread limits.

At exactly run396 seed9261703's best vector, original loss is
0.1357233047033634 and fixed spectrum-branch loss is 0.1464466094067266.
Gradient norms are 1.77612e-8 and 1.91640e-8; least Hessian eigenvalues are
-1.4609760402e-9 and -1.5584236709e-9. Neither meets the frozen -1e-6 trigger.
Seventeen Hungarian assignments give a distinct-root gap of
0.2133883472436626. Duplicate slots carrying one eigenvalue are not distinct
branches. Strict gap makes the selected branch locally active; branch descent
also implies lower-envelope descent without requiring persistent assignment.

Independent read-only audit398 reconstructs the complete base native/compiled
lists, objectives and unitarity, checks saved feasible alternatives/gap
arithmetic, Hessian symmetry and all eigenpairs, orthonormal COLUMN eigenvectors
and signs. It does not independently recompute Hessian entries or solve the
assignment-gap problems. The saved compiled candidate has 12 CNOTs but maximum
off-diagonal entry 0.3535533934 and is invalid.

Independent FD399 executed once but failed its hash guard before any matrix:
a scope-wording edit changed its point-audit dependency after manifest freeze.
It is closed failed with zero evaluations; no hidden retry occurred. Historical
metadata was explicitly reconstructed into a separate copy with exactly the
expected SHA-256, not represented as an original archived snapshot.
Replacement400 checked all dependencies before its separate reservation and
completed ten independent NumPy matrices once in 0.0667 seconds, under a
30-second one-thread cap. It tests only the two frozen least-eigenvector
columns, matching objectives, and signed steps 0.01 and 0.001. No AD, fit,
reassignment, or full independent Hessian was computed.

The saved FD second differences are respectively 2.31016e-6 / 2.23710e-8 and
5.42858e-6 / 5.34295e-8. No signed sample decreased its objective. Roughly
100-fold shrink for a ten-fold smaller step is consistent with O(h^2)
finite-step truncation. Finite-step effects and cancellation limit resolution;
the tiny negative AD signs are not independently resolved. The immutable
ledger400 notes and board102 outcome earlier attributed this to cancellation
alone; this report and FD400_CLOSEOUT.md qualify that unsupported attribution.

Runs398/400 are terminal inconclusive. No trigger-based escape fit, exact
certificate, incumbent change, local-minimum proof, or family exclusion follows.

## Collaboration and next hypothesis

Three GPT-6 Luna collaborators supplied algebra, static numerical review, and
independent verification, with direct peer findings and persisted handoffs.
Boards92/99 cover timeout review;96 verification;97 branch curvature;100 the
Bell proposal. Board95 was mistakenly claimed and closed without work;94 is an
unexecuted duplicate proposal, not a computation.

The best next structural hypothesis moves opposite-pair Bell decoders before
the fan-in, retaining all98 free coordinates:

`L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.

For F=exp(i aXX+i bYY), the exact decoder prefix uses input control
Ry(pi/2)S-dagger and target Rx(-pi/2), F(-pi/4,0), then X on control.
Place both pairs' input locals in L0, X0 only in L1 and X1 only in L2.
This exposes the Bell-coordinate target (CZ_A tensor I_B)SWAP_AB for
A=(0,2),B=(1,3). Suffix synthesis remains unproved. The prefix is an
initialization, not a frozen constraint; any next fit must release all98.
The reordered word is a different family and has not yet been computed.

## Verification and limitations

Saved FD arithmetic/count/no-decrease assertions were independently checked
without new matrices. Ledger and board audits passed. Focused policy/checker
tests passed (16 tests); adversarial review required the qualified FD
interpretation above. All run inputs, outputs, plans, candidate and failed
attempts are retained in experiments/runs/397--400. Exact acceptance remains
required for any complete passing gate list. The continuing success goal is
active, with established interval 6 <= C_min <= 13.
