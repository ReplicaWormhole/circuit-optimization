# Numerical round3: run373 / hypothesis21

Two distinct12-CX topologies were frozen after direct algebraic exchange, replacing source slots0–2 or6–8 by identity,12,02. All168 local parameters remained free; canonical13 local warmstart plus seeded0.08 noise, seeds926301/926302. L-BFGS-B analytic Torch gradients;300iterations/120seconds each,240seconds total, one thread. Run373 reserved and linked before fitting. Reproduce with `python3 collaboration/work/numerical_round3/run.py`; immutable plan/source/dependency hashes are in plan.json.

Both saved full gate lists failed independently reconstructed NumPy diagonalization: case0 maxoff0.4487654196973814, case1 maxoff0.6728896343792273. Unitarity errors below1e-15;12CNOT each. Fits stopped at130/157iterations,5.864232516seconds total. CLI check_circuit.py independently rejected both (expected exit1), and verifier_round3 confirmed the same values. No passing candidate and no exact-certification job warranted. Run373 failed; board21 inconclusive. This does not exclude either topology or change6<=Cmin<=13.

## Correction to architectural motivation

The frozen plan and actual chats are preserved verbatim. The proposed precursor transformation of source[01,32,12] into decorated braid[01,12,01] needs TWO changes:32->12 and last12->01. Earlier messages claiming32->01 alone suffice are wrong; that produces[01,01,12]. Actual searched12 schedules[12,02] are unchanged, still well-defined distinct refits. No equality with the source was claimed or used.

Recommended next bounded hypothesis: first refit13-CX schedules containing true[01,12,01] motif, then seek commutant-constrained local decorations permitting exact reduction. No follow-up run dispatched.

Adversarial review: corrected the precursor error above. Input conventionMSB qubit0, chronological gates, unrestricted eigenbasis/order. Numerical failure is bounded evidence only. No dependencies or shared checker edits. All own JSON/scripts/chat notes retained as intentional evidence; coordinator owns commits.
