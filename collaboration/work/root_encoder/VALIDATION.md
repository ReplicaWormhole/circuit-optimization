# Constructive round10 validation and review

Artifacts: math, computation, research scripts and workflow. Presearch verifier
binds exactP source, catalog and conservative costs. Beam381's integer action
is independently checked at every input; it is exactP, not a diagonalizer.
Compiler focused Fraction phase identities and closed-parity proof reviewed
independently. New full native67CX circuit is separately certified by existing
independent exact_check in run382, conductor64, correctUV4=DU and multiplicities.
Numerical check independently passes. Source and compiler hashes are retained;
old381 catalog and all prior evidence remain unchanged.

Adversarial review used `~/.codex/adversarial_review.md`. Corrected an overly
restrictive next-round statement: constant-on-orbit accumulated phases suffice
to preserve the current intertwining, but are not necessary when new eigenlabel
order is allowed. Future substitutions need explicit phase-aware checking.
No remaining defect in sign, controlled branch/global phase, Fourier bit reversal, overwritten-wire
encoder semantics, source binding or cost interpretation. The complete gate
list is explicitly larger than incumbent13; generic/family improvement is
not a smaller global upper bound or proof of optimum. Beam success gives an
upper bound, no minimum. Globalphase omission is common to every input, not
a silently discarded conditional branch scalar. No incumbent promotion,
shared checker changes, new dependencies, duplicate certificate or remote work.

New focused checks are exact16basis action, phasepolynomial identities, native
gate numerical and exact certificates, script compilation, hash/database and
diff audits. Unchanged shared22-test infrastructure suite need not be repeated.
Final ledger audit382runs196candidates no problems, zero running. Board audit
four submission hashes/JSON valid, not proof certification. Both SQLite
integrity checks return ok. Certificate source/compiler hashes match and
new algebraic/search/verifier scripts pass py_compile. Working and staged
whitespace checks pass. All workers complete, no scientific job scheduled.
Only assigned files/evidence included. Research JSON/native gates/scripts/notes/
chats and code/candidate/database snapshots are kept evidence; caches ignored.
Goalactive: explicit optimum and matching lower-bound evidence still missing.
