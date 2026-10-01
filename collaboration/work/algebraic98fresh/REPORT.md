# Fresh full98 initialization: structural scope

Claimed board hypothesis 73. This is a read-only structural handoff; no
numerical or symbolic search, fit, certificate, or candidate edit was made.

## Coordinate and family scope

The frozen run389 chart specifies 98 real coordinates: 84 SU(2) rotation-vector
coordinates for local layers `[L0prime,L1prime,L2,L3,L4,L5,L6]` (three
coordinates per qubit, q=0..3 inside each layer), six rotation-vector
coordinates for the retained target-wire gates `A1(q1), B1(q1)`, then eight
F-angle coordinates for pairs `(0,2),(1,3),(0,1),(2,3)`. Each local vector
`(x,y,z)` represents `exp[-i(xX+yY+zZ)/2]`; each F block is
`exp[i a XX+i b YY]`. The gate chronology is

`L0prime CX01 A1 CX21 B1 CX31 L1prime F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`.

This is the same relaxed family derived in `algebraic98/REPORT.md`: the
control-wire pieces of the original A and B layers are absorbed into the free
boundary layers, while A1 and B1 remain interior because q1 is the common CX
target. It retains four native CX and four F blocks, which compile to 12 CX
with arbitrary one-qubit gates. It is a restricted architecture family, not a
claim about every possible 12-CX circuit.

## Distinguishing the fresh point and interpreting it

Runs385/389 use the embedding with `A1=B1=I`, so their six A1/B1 rotation
coordinates are zero. Requiring both serialized A1 and B1 to be nonidentity
in the frozen initial vector guarantees a different coordinate point, and
therefore exercises interior target-wire freedom omitted by that embedding.
Freeze the exact 98-vector, parameter order/chart, initialization seed and
sampling rule before run reservation. Coordinate difference alone does not
show a distinct circuit modulo gauge/redundancy, a separate attraction basin,
or a likely successful fit; it only prevents repeating the inherited start.

A fit that reaches a verified diagonalizer would be a candidate in this family,
with native cost 4 CX+4 F and compiled cost 12 CX; it would still need
independent gate-list validation and an exact certificate before acceptance.
A failed bounded fit is optimizer evidence at one start, not exclusion of the
family, not a topology lower bound, and not a proof against 12 CX. Run389's
near-zero least curvature `-2.70e-9` did not meet its `-1e-6` trigger, while
57 Hessian eigenvalues were within `1e-6` of zero; it establishes neither a
local minimum nor a saddle and does not justify a same-point repeat.

## Handoffs

The numerical collaborator should freeze a reproducible full98 point with
explicitly nonzero A1 and B1, preserve all 98 coordinates, and keep the
single-start budget and interpretation bounded. The independent verifier
should check the exact chart/order, nonidentity A1/B1 in the serialized start,
chronology, and native-to-compiled matrix equivalence before or after the fit;
any candidate still requires its separate exact-certification process.
