# Frozen plan: 16 literal phase-aware variants

Board 117. This directory contains only the prepared screen; no circuit matrix
or objective has been evaluated. Root owns ledger reservation, run-workspace
creation, and authorization to execute.

The script enumerates the Cartesian product of four binary choices: P RCCX
control order and adjoint, and final negative-control CCH RCCX control order
and adjoint. Each candidate has the exact Bell decoder (2 CX), difference
router (2 CX), `CZ_A` replacement for `CS†_A` (1 CX), the P RCCX (3 CX), exact
controlled-H (1 CX), and target-conjugated CCH RCCX (3 CX), totaling 12 CX.
The script writes every complete gate list before evaluating it, then records
the normalized squared Frobenius off-diagonal loss, maximum off-diagonal
entry, unitarity/root errors, CX count, and content hash. A numerically passing
best candidate is retained and sent once to `exact_check.py`; no candidate is
promoted automatically.

The relative-phase labels do not imply validity. The screen checks each full
unitary against the right-shift diagonalization condition. A failed finite
screen is bounded evidence only and does not exclude the 12-CX family.

Frozen limits: exactly 16 variants, deterministic order, one thread, no
optimizer, no random draws, 30-second command limit, and one exact check only
if the numerical checker accepts a candidate. The run wrapper must reserve an
experiment ID and copy this directory into `experiments/runs/<ID>/` before
execution.
