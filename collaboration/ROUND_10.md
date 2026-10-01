# Exact encoder and a compiled orbit construction

Initial clean d44ffdb; live owned-agent registry confirmed three completed
workers, process inspection no scientific job, ledger380 zero running, board
43–45 terminal. Prior turn is progress: tested architecture hypotheses and
established a concrete orbit family for synthesis. Full optimum goal remains
unproved with6 <= C_min <= 13. New board46–50 and isolated encoder directories
preserve all old evidence. No incumbent certificate work is duplicated.

## Exact encoder search

Algebraic board46 derives all four Boolean ANFs of fixed orbit-address P
and a44-CCCX reversible baseline. It checks overwritten-wire semantics by
the complete16row action, not by assuming simultaneous ANFs are a circuit.
Independent verifier checks ANFs, rows and catalog before any beam search.

Numerical board48/run381 freezes exact P, source report hash and108 mixed-
polarity X/CX/CCX/CCCX catalog. CompiledCX upper weights0/1/6/34, local gates
and native depth separate; repeated-state checks prevent free-X loops.
Deterministic beam128, depth12, expansion200000, wall120s, one thread.
It finds EXACT P in0.637379230s, depth10,122364expansions:

1. CX3→2
2. CX0→1
3. CX2→1
4. CX2→3
5. CX3→0
6. CX1→0
7. CCX controls(q0=1,q3=0), targetq2
8. CX1→0
9. CCX controls(q0=1,q2=0), targetq3
10. CCX controls(q1=1,q3=1), targetq0

Independent integer action verifies all16basis rows, q0 MSB and chronological
gates: `[0,4,7,8,6,2,11,12,5,9,3,13,10,14,15,1]`. SevenCX and threeCCX give
compiledCX upper25. Negative controls can use four free local-X conjugations;
the later phase compiler handles their signs directly. Finite beam is not
a proof of minimum encoder cost. Encoder alone does not diagonalize V4.

## Explicit full gate list and independent exact certificate

Root board49 derives shared Gray-code parity compilation for an m-bit projector
phase at2^m-2CX, instead of compiling each Pauli parity independently. Gray
cycles close each target, use no ancilla, and omit only an input-independent
global scalar. Assigned verifier implements it and performs focused exact
Fraction phase-polynomial checks. This givesCCX6,CCCX14 and conditionalH
cost6/14 for two/three controls. Old run381's frozenCCCX34 catalog is not
rewritten; improved catalog costs require a future reservation.

Global bit-reversedF4 on the two low bits costs2CX. Conditional inverse of
the SAME transform in the00upper sector costs26; three-controlledH for the
pair sector costs14. The mixing block costs42, not the previous loose268.
Combined with encoder25, exported complete source is
`work/verifier_label12/encoder/orbit_encoder67.json`,67CX and rational-pi local
gates. Negative controls, conditioned phases and Fourier orientation are
explicitly implemented; there are no free native multi-control primitives.

Root board50/run382 reserves one independent `exact_check.py` invocation,
180seconds, one thread, source/compiler hashes frozen. It returns exit0 in
under one second: exact_diagonalizer true, cyclotomic conductor64,67CX, roots
counts(6,3,4,3). Labels are
`[0,0,0,2,0,2,1,3,0,2,1,3,0,2,1,3]`, roots i^label. This certifies the full
normalized native gate-list UV4=DU exactly, separately from the encoder's
truth table and the compiler's phase identity. Numerical checker agrees,
maxoff5.79e-16, unitarity1.55e-15. Certificate output/source binding is in
`work/root_encoder/`. No exact13 certificate was rerun or promotion made.

The67CX construction is NOT a reduced incumbent: accepted exact13 remains
better and the global lower bound remains6. It supplies a concrete jointly
optimizable construction in place of the3838CX abstract compilation bound.
It is not total-spin diagonalization or strong Schur. Separate nativeCX/F
model retains6CX+4 independently tunableXX/YY gates, depth8, compiledCX upper14
and numerical native JSON with analytic provenance; no10CX implication.

Recommended next round, not dispatched: exact phase-aware joint simplification
of this native gate list, rather than only permutation synthesis. First apply
local/CX commuting and adjacent-inverse cancellations, preserving the complete
matrix and recording gates/counts. A relative-phase CCX substitution requires
explicit accumulated phases: constant-on-cycle-orbit input phases suffice to
preserve UV=DU. Otherwise recheck diagonalization with arbitrary new label
order, or supply explicit compensation; constancy is not a necessary condition.
Any candidate must be native,
independently checked and exactly certified; a count above13 is family-local
improvement only. Bound any enumeration to one frozen rewrite pass/template
and a stated resource limit before reservation; do not silently repeat beam.

Actual direct chats and board handoffs are in `CHATS.md` and each researcher
encoder directory. Scientific ledger382 and all hypotheses are terminal;
new scripts/JSON/reports/databases are intentional retained research evidence.
Final source/integrity/compiler/diff audit is in `work/root_encoder/VALIDATION.md`.
