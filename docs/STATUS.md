# Current research status

Reviewed on 2026-10-01 through run421 and stopped round24 preparation. Check the live
[scientific ledger](../experiments.sqlite3) and
[board](../collaboration/board.sqlite3) before starting new work; this page is
a reviewed summary, not a process monitor.

## Established result

In the ancilla-free, all-to-all model with arbitrary one-qubit gates and
directed CNOTs, **6 <= C_min <= 11**. The
[11-CNOT gate list](../experiments/runs/411/candidate.json) has 30 elementary
gates and exact rational-pi angles, with an
[acceptance record](../collaboration/EXACT11_ACCEPTANCE.json). Two independently
implemented exact conductor32 certificates pass; run412 explicitly checks all
256 diagonalization and 256 unitarity identities. The exact12, exact13 and
exact14 constructions and acceptance records are preserved. Optimality remains
unknown; total-spin diagonalization and a strong Schur transform are unclaimed.

## Latest reviewed work

Work is paused at the user's request. Round24 produced a by-hand rational-pi
24-primitive/four-CX solution for the custom reduced v=0 target CZ02 CX10,
independently reviewed for complete-word phases. Computational exact certification
is pending; this is not a full V4 diagonalizer. The obsolete reduced numerical
draft and proposed four-axis full10CX extension remain unexecuted. All owned
boards are terminal; no new ledger attempt beyond421 or active scientific process.
See [the stopped-round handoff](../collaboration/HANDOFF.md) and
[round24 preparation](../collaboration/ROUND_24.md). Resume only on explicit
instruction with fresh claims, frozen code/inputs/budgets, and reservations.

Round23 tested the reordered ten-CX schedule
`02,13,01,23 | 02,32,10,02,12,32`, with all 132 local coordinates free.
Preflight 418 stopped before unitary/RNG/AD evaluation on a descriptive-metadata
comparison; corrected 419 passed 22 full unitary builds, 3 AD and 6 FD checks.
One bounded two-start fit, run 420, completed 136 iterations and 152 calls in 0.942 seconds.
Best by loss is case0: loss 0.125, maxoff 0.7071; case1 has loss 0.26278,
maxoff 0.6361. Both are invalid. Independent endpoint 421 passed 14 full matrices
and complete word/vector/noise/topology/history/metric checks. All runs closed;
exact11 retained and no topology or local-minimum exclusion follows. See
[round23](../collaboration/ROUND_23.md). The elapsed-time error in421's ledger
notes is documented separately; saved result/stdout record 0.173786 seconds.

The [new restricted phase argument](../collaboration/work/algebraic10reorder/PHASE_OBSTRUCTION.md)
shows why a literal early-CH transplant cannot work with computational control
bases fixed and no later x-basis mixing. The full 132-coordinate family allows those changes
and remains open. Best next hypothesis: solve the exact v=0 operators W00=CZ
and W10=[[0,I],[Z,0]] jointly with tilted x controls and a later x readout,
then extend to v=1. No phase-compatible angles or ten-CNOT word are yet known;
no next computation is dispatched without a fresh finite freeze.

## Previous ten-CNOT deletion investigation

Round22 tested ten-CNOT routes from the certified11 circuit. Reservations413
and414 stopped before scientific execution for board ownership and ledger-target
guard errors; both are closed and preserved. Corrected independent preflight415
passed121 full matrices,11AD gradients and22finite-difference checks. Sole bounded
fit416 completed all11 single-CX deletions with132 free local coordinates,
834iterations/938calls in6.171seconds. Best deletion7 has loss0.013675822916672686
and maxoff0.2044593930001434, invalid1e-9. Independent endpoint417 passed77 full
matrices and all endpoint/word/noise/topology/history/hash checks. All runs closed;
no valid10 circuit, exactness check, incumbent change or topology exclusion.
See [round22](../collaboration/ROUND_22.md).

The [restricted suffix argument](../collaboration/work/algebraic10/SUFFIX_OBSTRUCTION.md)
rules out saving one target interaction with the stated prefix and sector labels
fixed. It leaves changed earlier blocks and general output-label freedom open.
The [next unexecuted hypothesis](../collaboration/work/algebraic10/NEXT_HYPOTHESIS.md)
moves CH earlier inside the merge: `02,13,01,23 | 02,32,10,02,12,32`, seeking a
six-CX suffix through joint wrapper/phase synthesis. A bare CNOT triangle permits
absorption, but the wrapped circuit identity remains unproved. No follow-up run
is dispatched; new computation requires a separately frozen finite budget.

## Previous eleven-CNOT milestone

Runs408–412 continued from the accepted12 circuit. Independent preflight408
passed all twelve deletion initializations and gradient checks. Bounded fit409
finished twelve fully free local11CX topologies; none met 1e-9, but deletion9
reached maxoff5.06e-8. Independent complete endpoint audit410 passed. The near-hit
inspired a by-hand construction that rotates the merged output basis and replaces
the three-CNOT selector with two controlled reflections, preserving all phases.
Complete-list numerical and exact411 passed; independent Fraction Phi32 run412
proves all 256 UV4-DU and 256 UUdag-I entries exactly zero. All attempts are closed;
Newton preparation was superseded without execution. See
[round21](../collaboration/ROUND_21.md) and the
[derivation](../collaboration/work/algebraic11nearhit/EXACT11_PROPOSAL.md).

Original fixed 98-coordinate word membership remains unresolved. Round22 above
tests the subsequent joint reflection/selector hypothesis within finite scope.

## Previous twelve-CNOT milestone

Runs405–407 developed and certified a Bell-label difference construction.
Run405 validates the primitive19CX comparator. The three GPT-6.1 Sol
collaborators repaired selector phases, independently checked a four-CX merge,
and froze the complete 39-gate, 12-CX candidate before computation. Run406's single
numerical diagnostic gives maxoff5.67e-16 and its exact cyclotomic check passes.
Independent run407 uses Fraction arithmetic modulo z^8+1, importing no existing
checker, NumPy or SymPy; all 256 UV4-DU and all 256 UUdag-I entries are zero.
All runs are closed. See [round20](../collaboration/ROUND_20.md).

This establishes an exact 12-CNOT diagonalizer. Membership of this circuit in the
original fixed 98-coordinate word remains open. The best next structural
hypothesis is to absorb a difference-router interaction into the merged
controlled block using output computational-permutation freedom, seeking 11 CNOTs.
No follow-up computation is dispatched; any search needs a new finite budget.

## Earlier bounded full98 investigations

Run398 completed both full98 Hessians at run396 seed9261703's repeated
plateau. Neither least eigenvalue (about -1.5e-9) meets the frozen -1e-6
escape trigger. Independent gate-list/eigenpair audit passed; run400's ten
independent directional finite-difference samples found no decrease but do
not resolve the tiny AD signs. Finite-step effects and cancellation limit
resolution. Run397 was stopped before execution for a timeout-handling defect;
run399 failed a stale source hash before matrix evaluation. All four attempts
are closed and preserved. See [round18](../collaboration/ROUND_18.md).
No local minimum, family exclusion, exact certificate or incumbent change follows.

Run402 tested the reordered Bell-decoder full98 word from one frozen
Bell-derived start plus seed9261801 Gaussian perturbation0.1 on all98. All98
were free. Independent preflight404 passed three-point matrix/serializer and
prefix checks. The fit completed113iterations127calls in1.365seconds, best
loss0.375 and maxoff1, invalid compiled12CX. Run401 was unexecuted for a
missing config key; preflight403 failed a runtime-version guard before matrices.
All attempts are closed; independent full initial/best gate-list audit passed
and confirms invalidity. See
[round19](../collaboration/ROUND_19.md). No certificate or exclusion follows.

The next structural hypothesis constructs the Bell-label exchange eigenbasis
suffix initialization, since placing a valid decoder prefix does not solve the
remaining transform. Any numerical test needs a fresh frozen finite budget and
reservation. At that round the exact12-CNOT success goal was still open; round20 resolves circuit existence.

Runs395/396 previously tested the spectrum-assignment objective with all98
coordinates free. Eight fresh starts returned to labeled loss0.1464466094,
original loss0.1357233047, maxoff0.3535533934. Full endpoints independently
failed diagonalization. See [round17](../collaboration/ROUND_17.md).

Run394 completed 32 full98 starts across four local rotation scales in
117.314 seconds, within its frozen single-threaded budget. The best seed9261616
has normalized loss0.125 and maximum off-diagonal entry0.7071067812.
Independent NumPy/SciPy verification checked all64 endpoint native/compiled
gate lists and the globalbest lists; every compilation has12 CNOTs, but none
diagonalizes V4. The run is closed failed, with no exact certificate or incumbent
change. See [round16](../collaboration/ROUND_16.md).

An independently checked opposite-pair Bell-label derivation gives eigenvalue
multiplicities(6,4,3,3), but its synthesis in the frozen word remains open.
The proposed eigenvalue-assignment residual was tested in395/396 above.
The 12-CNOT goal was open at that round; round20 resolves it.

One fresh full98 initialization with nonzero interior target-wire gates was
fitted in run393: seed9261501, all98 coordinates free, one thread, maximum300
iterations and external120-second timeout, no restarts. After284 iterations,
the normalized loss fell from0.9085 to0.4549; the best compiled12-CNOT circuit
still fails diagonalization (maximum off-diagonal entry0.67436). Independent
NumPy/SciPy checks validate all initial/best native and compiled gate lists.
No exact certificate or incumbent change. See [round15](../collaboration/ROUND_15.md).

The earlier full98 Hessian/independent derivative checks failed the frozen
negative-curvature trigger; see [round14](../collaboration/ROUND_14.md).
These bounded results prove neither a local minimum nor family exclusion.
No fit is currently dispatched. A new seed or other hypothesis requires a new
frozen finite budget and ledger reservation.

## Where to inspect evidence

- [Accepted and baseline circuits](../collaboration/incumbents.json)
- [Round reports](../collaboration/README.md) and [actual recorded exchanges](../collaboration/CHATS.md)
- [Historical search chronology](history/SEARCH_14_STATUS.md)
- [Handoff notes](../collaboration/HANDOFF.md)

Use `python3 experiment_log.py summary` and `python3 collaboration_board.py list`
for live local status. A numerical fit or a board submission is not an exact
certificate or a global lower bound.
