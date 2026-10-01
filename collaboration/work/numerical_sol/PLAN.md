# Frozen literal 12-CX phase screen

Board 120; numerical_sol owns preparation, coordinator owns reservation,
execution, finish and commit. Existing numerical98phase12screen artifacts are
preserved. This task touches code, math, computation and workflow.

The exact chronological Bell prefix is CX02,H0,CX13,H1 (2 CX), followed by
R=CX01,CX23 (2 CX), then CZ_A=H2,CX02,H2 (1 CX). P is the literal 3-CX RCCX
on controls (3,0), target2. Exact CH(1->0) is Ry0(-pi/4),H0,CX10,H0,Ry0(pi/4).
The final relative-phase CCH is X1,Ry2(pi/4),S2,RCCX(1,3;2),Sdg2,Ry2(-pi/4),X1
(3 CX). Each RCCX independently chooses either control order and either the
forward word or its literal reversed inverse. Full RCCX chronology is
H_t,T_t,CX_b,t,Tdg_t,CX_a,t,T_t,CX_b,t,Tdg_t,H_t; T/Tdg use phase-only U3.
This yields exactly16 declared complete words and12 CX each. Literal adjoint
words differ, but may represent the same unitary; no labels are pruned.

Seed0, no random draws, complete NumPy unitary diagnostic, exactly16 cases,
one thread,30-second total in-process evaluation budget and coordinator
external30-second command cap. There is no optimizer and no exact certificate
inside this screen. The commander must reserve the exact config.json object
with method "complete NumPy unitary diagnostic", then link run to board120.
The script checks the matching active ledger config/code/method and board link
before matrices. Exclusive execution_started.json forbids silent retries.

Preparation command (serialization only):
`python3 collaboration/work/numerical_sol/screen.py --prepare`.
Coordinator copies screen.py,config.json,variants.json,PLAN.md to the new
experiments/runs/ID workspace. Execution command:
`timeout 30s python3 experiments/runs/ID/screen.py`.
All16 gate lists are saved before the first evaluation. Every case records
shared-checker errors, normalized off-diagonal loss, direct UV4-DU residual,
nearest diagonal-root multiplicities in order(1,-1,i,-i), hashes and equivalent
unitary representatives. Expected valid multiplicities are(6,4,3,3); nearest
counts for invalid cases are diagnostic only. Best complete candidate and
every result are retained; failed and incomplete attempts must be closed.

Qubit0 is most significant; all gate lists are chronological. No valid
diagonalizer is assumed. Algebraic_sol reported a by-hand obstruction in
sector d=(q1,q3)=(1,0),a1=1 for all frozen CZ variants, and RCCX Hermiticity
makes forward/adjoint choices redundant. The finite numerical screen may
corroborate that restricted claim; it does not exclude other12CX constructions,
certify a circuit exactly, or change the accepted13CX incumbent.
