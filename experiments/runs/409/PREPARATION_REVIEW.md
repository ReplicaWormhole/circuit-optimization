# Preparation review, no numerical execution

Artifacts touched: code, computation and workflow. Main script and family
were syntax-checked with `python3 -m py_compile` (pass). Pure JSON topology
generation verified12 sourceCX positions and12 schedules with11CX each.
Runtime/library versions and module hashes were collected without importing
the family, constructing a circuit matrix, drawing perturbations or fitting.

Static adversarial review found and fixed a timeout issue before freeze:
the total110-second alarm now uses a distinct BaseException class to bypass
case-level error handlers and preserve best vectors/words before stopping.
Case8-second stops remain ordinary bounded case outcomes. Saved bests do not
require a new matrix evaluation during terminal-exception cleanup.

The all144 free-coordinate mapping, chronological source-gate collapse,
determinant/global-phase handling, U3 singular chart, smooth sinc rotation
formula, CPU autograd, full12case seed allocation, all initial-word saving,
one-shot marker, exact config membership/value/hash guards and source topology
rule were inspected. No further confirmed defect found in static review.
Runtime matrix/serializer/gradient correctness remains unverified until the
separately reserved independent preflight passes. Any later change requires
new hashes and invalidates the current freeze.

The root and verifier received the frozen code/config hashes. The coordinator
owns reservation, preflight dispatch, fit dispatch, terminal ledger outcomes
and commits. No shared checker, earlier candidate or concurrent agent files
were modified. This new directory and JSON packet are retained evidence;
Python cache is ignored. No optimizer execution or result claim is made.

Algebraic_sol handoffs identify two structural alternatives not silently
tested here: CX01 can commute adjacent to CH, and conjugating the merged block
by routerCX23 exposes a phase-weighted controlled swap. Their costs are not
proved, and this run retains the original exact12 CX order after each deletion.
These alternatives therefore need separate hypotheses and budgets if pursued.
