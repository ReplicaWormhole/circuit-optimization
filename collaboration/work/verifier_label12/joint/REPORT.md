# Independent joint-label review

Hypothesis30; computation, math and workflow. No optimization or exact13
recertification. Full U3/CX matrices use the independent bit-column checker
from the preceding round, with q0 MSB and chronological UV=DU conventions.
The source unchanged hash is9a9dd0e0b430e935f99f4d998fb25d5e1b48160ad16a6c46a27e2574c59c2398.

All prescribed labels preserve root counts(6,3,4,3). The full source sequence
with identity entangler slot1 and all56 local gates is parsed sequentially;
optimizer reconstruction differs from the independent matrix by at most
1.4488835837920476e-15 up to global phase. Frozen target swaps and hashes
are independently verified in independent_audit.json.

First composed swap(2,6)(3,7):12CX, maxoff0.1517277915509333,
offdiagonal squared loss/16=0.00942310839202986; selected-target loss/16
0.009459905009667346, maxUV-DU0.0689367816905082. This is an improved invalid
lead relative to source maxoff0.2588190369994409/loss0.050240473580837317.
Second swap(8,13)(9,12):12CX, maxoff0.7071067918307067, offdiagonal loss0.375,
selected-target loss0.46343390751450736. Both lists are invalid, with
unitarity errors below1.12e-15. No exact incumbent promotion follows.

Independent Gaussian-integer sums reproduce reduced q3 spectra: base and
firstjoint{0,2}, secondjoint{-2-2i,4+2i}, each secondsingle{-1-i,3+i}.
Secondjoint reduced q1 spectrum returns{1+i,1-i}; firstjoint gives{3-i,-1+i}.
Product-local output conjugation preserves reduced-target eigenvalues, so
the two joint targets are distinct under this equivalence and secondjoint
is distinct from its constituent singles. Its bit permutation toggles q1
and q3 iff q0=1 and q2=0. These exact structural statements concern target
labels, not circuit implementation costs or universal optimality.

The coordinator's halfshift-support consequence was separately reviewed:
with U V U†=D and A=U†p_jU, translating the existing theorem uses U_old=U†.
If A commutes with V², its minimal Pauli support is invariant under(02)(13).
The only invariant subsets are empty, the two opposite pairs, and full.
Tracelessness forbids empty; the existing opposite-pair theorem forbids both
pairs. Hence support is full. This imports completed exact evidence and
requires genuine Pauli-axis images of diagonalizers; no certificate rerun or
new exclusion of the3024 surviving six-CX schedules is claimed.

Review of frozen budget: two200iteration L-BFGS fits and at most5SVD steps
each,120seconds total, one thread; observed8.639925seconds. Guard checks
between operations do not preempt a single long Jacobian/SVD. Costs remain
fixedCX+arbitrarylocals: accepted13CX upper bound, candidate12CX counts
supply no upper bound while invalid. Separate hybridCX/F baseline10native
entanglers with compiledCX upper14 remains a different gate model.

Commands: verify.py with both saved targets writes independent_checks.json;
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3
collaboration/work/verifier_label12/joint/audit.py passes. Adversarial review
against ~/.codex/adversarial_review.md found no confirmed defect. Hashes,
chronology, phase freedom, selected versus nearest targets and bounded-failure
scope are retained. No exact claim, global optimum or topology exclusion.
All scripts/reports/JSON/chats are intentional retained evidence; coordinator
owns shared files and local commit. No push.
