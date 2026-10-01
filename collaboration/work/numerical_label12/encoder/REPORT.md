# Exact reversible orbit encoder from bounded beam

Board48, ledger381. Source orbit-address report hash and exact permutation
table frozen before search. Assigned verifier independently approved all108
catalog permutations/cost bounds and table before reservation. This is an
encoder-only combinatorial search, not a whole-circuit diagonalizer fit.

Beamwidth128, depthcap12, expansioncap200000, wallcap120seconds,1thread,
no randomness. Catalog mixed-polarity X/CX/CCX/CCCX with CX compiled upper
weights0/1/6/34. Local X count separate. X has unit native depth and seen-state
cost/depth comparison prevents free-gate loops. Rank and deterministic ties
fully stated in config; score is a heuristic, not an optimality proof.

Found EXACT P after10native gates,122364expansions,0.637379230seconds. Rows:
`[0,4,7,8,6,2,11,12,5,9,3,13,10,14,15,1]`, q0 most significant. Native
chronological steps:

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

SevenCX plus threeCCX compile to at most25CX via exact standard six-CX
Toffoli decomposition. Two negative controls require four free local-X
conjugations; noCCCX occurs. Count25 is an encoder upper bound, not a minimum
and not a new V4 diagonalizer. Independent verifier checks all16basis rows,
chronology, gate action and cost separately. Fullnative list and trace retained
in result; no approximate encoder is promoted. No standard gate-list JSON was
submitted to cycle checker because it does not describe U diagonalizing V4.

Adversarial review: checked rightcycle table orientation, bijection, catalog
allpolarities/bitorder, chronologicalcomposition, deterministic scoring/ties,
actualresource caps and separate native/CXupper/local costs. No confirmed
defect. Finitebeam/dedup can miss betterpaths; exactP does not imply optimality.
CX25 alone exceeds13 incumbent, and conditionalFourier cost remains additional.
Keep result as constructive encoder improvement over generic44CCCX1496bound.
Next distinct step is analytic joint encoder/Fourier synthesis with explicit
compiled gates and cancellations, under a new budget; no further search run.
Research files kept, Python caches ignored, coordinator owns commit.

Independent verifier has now confirmed all16basis rows and exact catalog
cost recount in `verifier_label12/encoder/beam_check.json`. Board48 completed.
Frozen CCCX34 catalog remains unchanged; a subsequently proved better CCCX
compiler bound belongs to a future catalog. No additional beam authorized.
