# Exact12-source deletions and router absorption

The previous turn made concrete progress: runs406/407 independently certified
an exact12-CNOT diagonalizer, accepted in EXACT12_ACCEPTANCE.json and committed
as14f7e36. This continuation seeks a smaller circuit without changing that
incumbent unless a complete independently validated gate list has an exact
certificate. Original full98-word membership is a separate open question.

## Frozen hypothesis before computation

The structural hypothesis is to absorb a difference-register router CNOT into
the adjacent controlled-reflection block, allowing arbitrary computational
output permutations and bases in degenerate eigenspaces. Algebraic11router and
verifier11router have distinct by-hand assignments and exchange findings with
the numerical collaborator directly and on the board.

The numerical diagnostic starts from the newly accepted12 circuit, not a prior
13-CNOT source. Delete each CNOT in chronological order0..11, leaving11 CNOTs,
and collapse intervening local gates into12 layers of four arbitrary one-qubit
rotations (144 free real coordinates). Initialize from the exact source after
deletion and perturb every coordinate by Gaussian sigma0.025radians using
NumPy PCG64 seed10012101 in fixed deletion order. One L-BFGS-B fit per deletion,
at most80 iterations each, one thread, total internal110s/external120s. No
restarts or deepening. Every base/initial/best complete gate list and parameter
vector must be retained. Source, dependencies, method, seed, topology rule and
limits are frozen in the ledger before execution; independent preflight checks
must pass first.

## Outcome

Run408 independently checked all twelve initializations:132 full matrix builds,
12 directional autograd gradients and24 finite-difference losses. It passed
before fitting (max phase error1.78e-15, gradient discrepancy1.34e-9).

Run409 completed all12 fits in7.108s,951 iterations/1067 evaluations. Best
deletion9 (first CX32 in the final selector) reaches loss2.27947e-15 and
maxoff5.05805e-8; deletion7 also reaches maxoff8.98710e-8. Neither passes
the1e-9 validity tolerance. The run is closed failed/numerical, preserving all
endpoints. Independent run410 checked84 complete matrices, all base/initial/
best vectors, gate lists, histories, topology and seed allocations; it passed,
confirming the invalid near-hits without a topology exclusion or local-minimum
claim. A separately budgeted Newton preparation was superseded before execution.

## New exact eleven-CNOT hypothesis before certification

The near-hit inspired a by-hand phase-compatible construction. Change the last
merged reflection from N4=(-X+Z)/sqrt(2) to
N4'=Ry(-pi/4)N4=-cos(pi/8)X+sin(pi/8)Z, keeping four CNOTs in the merge.
This adds a controlled-v Ry_y(-pi/4) to the old merged unitary. After controlled-H,
the v=1 sectors have target axes H=(X+Z)/sqrt(2) and M=(Z-X)/sqrt(2).
Chronological CZ(u,y), then controlled-v Nv, Nv=sin(pi/8)X+cos(pi/8)Z,
maps both to Z, replacing the old three-CNOT final selector with two.

The independent verifier checked all four sectors including scalar phases:
CZ_xy; (-i)^y Z_x; i^x Z_y; f_i(y)Z_x, with f_i(0)=i,f_i(1)=1.
Predicted labels in roots(1,i,-1,-i):
[0,0,0,2,0,1,3,0,0,1,2,3,2,3,1,2].
The full circuit cost is2+2+4+1+2=11, with exact rational-pi angles and no
ancilla. Two independently implemented complete gate-list exact certificates
over conductor32 were frozen before computation and passed as described below.
The exact12 construction and acceptance record remain preserved.

## Certification and acceptance

Run411 passed its single numerical60s and exact90s checker calls in0.0752s and
0.4467s, one thread, seed0/no randomness, external160s cap. Complete30gate11CX
word, depth8, numerical maxoff5.03934e-16; exactPhi32 labels match the prediction.
Frozen scriptSHA dc3686e6ee071789c489c7b8f419c2fcbd65a6e5e97f7642b7445e232e1315ea;
candidateSHA557c55cddefb697bcd4e038db300117d0aaaf6f6af97e851e0831baf257dbe4a.
Its entire gate word equals the algebraic collaborator's independent serialization.

Run412 passed independent Fraction arithmetic in Q[z]/(z^16+1),z=exp(i*pi/16),
all256 UV4-DU and256 UUdag-I entries zero. One90s invocation, one thread,
seed0/no randomness,0.1764s; discovered labels match with multiplicities6/3/4/3
in roots1/i/-1/-i. SourceSHA
f0c9c05f872139937e1f9fd5ad22006033ec1d3f7e3467788f4d1199b4b6eb09;
saved full exact matrixSHA
fbb007935825221f250efc5d86fb0dd3db1e438f92706ef7eac0ae457a81e965.
No existing checker, NumPy or SymPy import. Both runs closed complete/symbolic.
Each workspace retains inputs, config/dependency hashes, command, outputs,
complete sole/best candidate and scope. Standard ledger ingestion/audit also
recomputes floating metrics as bookkeeping validation, without search.

Independent source/full-list/output/hash review passed in
work/verifier11router/ACCEPTANCE_REVIEW11.md and work/numerical11exact/REPORT.md.
Acceptance is EXACT11_ACCEPTANCE.json. Established interval6<=Cmin<=11; the
preserved6-CX lower bound is inherited, not re-proved. No optimality, novelty,
full98-membership, total-spin or strong-Schur claim.

## Next hypothesis and closeout

Next seek10CX by combining the last merged reflection with the two final
selectors, allowing conditional target-basis changes. This is untested and no
next computation is dispatched. A new diagnostic needs its own finite budget,
board ownership, frozen inputs/code/method/seed/dependencies and ledger row.

Adversarial review found no remaining blocker in the actual gate list or exact
certificates. A misleading sentence about Ry's axis in the derivation was
corrected; equations and candidate bytes were unaffected. Twenty focused board,
ledger and checker tests passed; existing SQLite ResourceWarnings did not fail
tests. Ledger/board audits passed. Historical rows/events from HEAD14f7e36 and
exact12/13/14 source certificates remained unchanged. All new scripts, JSON,
matrices, logs, reports, databases and unexecuted Newton preparation are retained
research evidence. Caches/journals ignored; no dependencies or remote operations.
Incumbent index, README, STATUS, artifact index and local policy identify exact11.
