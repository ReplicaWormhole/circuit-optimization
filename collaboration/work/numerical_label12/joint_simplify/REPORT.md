# Exact-preserving local rewrite of orbit circuit

Board51/run383, exact parent382 source orbit_encoder67.json SHA
d595df83cc218ab4a532cdc2c51ab90a19dd68503fccc8d10311a23f0608ed6d.
Assigned verifier confirmed source and frozen rewrite script before reserve.
One bounded deterministic first chronological i/j rewrite,100sweeps,
50000rulechecks,120seconds,1thread. No optimizer, beam, relative-phase
Toffoli replacement, numerical U3 merge or additional rules.

Rules: disjoint supports commute; CX pairs commute iff neither control is
the other target; controlRz and targetRx commute with CX; same-axis same-wire
rotations commute. Identical H/CX cancel through these commutants; exact
Fraction rational-π same-axis rotations add without angle modulo reduction
or discarded scalar phase. Restart after each rewrite. Code/config/source
hashes archived in ledger before execution.

Fixpoint after1225rulechecks,0.000889951seconds:149→143total gates,67CX
unchanged, local count82→76. Five exact Rz merges: one zero pair removes two
gates, four nonzero additions each remove one. Complete trace and native
candidate saved; SHA99742270d05745155e74379602293bfb7bf2c19782d57def349521a570a000f5.
Numerical gate-list checker passes maxoff7.627804406540184e-16. Independent
verifier replays exact trace and compares matrices; root must separately
reserve full exact certificate of the changed JSON. No duplicate exact13 work.

This reduces local gate count but supplies no CX improvement over67 and no
incumbent improvement over13. It is not a global rewrite optimum or exclusion
of more powerful joint synthesis. Relative-phase encoder truth tables alone
are insufficient for diagonalization without reviewed orbit-phase constraints.

Adversarial review: source binding, gate chronology, commuting-crossing proof,
Fraction addition signs, no2πphase deletion, branch/cap limits and full trace
retained. No confirmed defect. Numerical pass is not exact certification;
await root's independent full exact check. All research evidence kept, cache
ignored, no shared edits, coordinator owns commit. No further computation.

Independent verifier subsequently replayed all five trace rules exactly and
checked full matrix agreement8.40e-16, preserving the candidate hash. Board51
completed; full exact UV certification remains assigned to root separately.
