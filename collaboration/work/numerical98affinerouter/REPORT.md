# Affine obstruction and conditional routing for Bell-label exchange

Board 113. This is a by-hand structural derivation only. I inspected the
existing Bell-prefix report and saved run402 outcome, but performed no matrix
evaluation, symbolic search, optimizer, or fit.

## Why a fixed affine router cannot use one orientation bit

Write a computational Bell-label pair as x=(a,b), where each label is a
two-bit string. Exchanging labels sends x=(a,b) to x'=(b,a). Their XOR
difference is

    Delta_d = x xor x' = (d,d),    d = a xor b.

The six distinct unordered label pairs realize the three nonzero differences
d=01,10,11, each twice. Any fixed sequence of CNOT and X gates is an
invertible affine map over GF(2), x -> Mx+c, with M invertible. Translation
c cancels in a pair difference, so the image difference is M Delta_d.

The three input differences are distinct nonzero vectors and
Delta_11 = Delta_01 xor Delta_10. Invertibility preserves distinctness, so
their images cannot all equal the same fixed unit bit vector. A stronger
obstruction holds: they cannot all have Hamming weight one even if the
orientation bit were allowed to vary by unordered pair. If
u=M Delta_01 and v=M Delta_10 were distinct unit vectors, then
M Delta_11 = u xor v has weight two. Thus a single affine CNOT/X permutation
cannot send every pair's two orientations to adjacent codes. This excludes
only that fixed-affine routing scheme.

## Conditional orientation routing after a difference-register transform

Use the interleaved wire order q0=a0, q1=b0, q2=a1, q3=b1. First apply
the two CNOTs q0 -> q1 and q2 -> q3. This difference transform D stores
a=(a0,a1) and d=(d0,d1)=a xor b in q1,q3. Swapping the original labels
leaves d unchanged and sends a -> a xor d.

For nonzero d, conditionally encode one unordered-pair coordinate t and one
orientation bit r as follows:

| Difference d | t | r | Effect of swapping labels |
|---|---|---|---|
| 01 | a0 | a1 | t unchanged; r toggles |
| 10 | a1 | a0 | t unchanged; r toggles |
| 11 | a0 xor a1 | a0 | t unchanged; r toggles |

For d=00, leave (a0,a1) as two singleton coordinates and do not mix them.
For each fixed d, the map from a to (t,r) is bijective. Together with d,
this defines a reversible conditional permutation. It is nonlinear as a
four-bit map because the choice of orientation coordinate depends on d, so it
escapes the affine obstruction.

Let sigma_x=(-1)^(x0*x1). The sign relation is constant in each routed
(d,t) block; the needed transform on r is:

| d | t | Sign relation | Transform on r |
|---|---:|---|---|
| 01 | 0 | equal (+,+) | H |
| 01 | 1 | opposite (+,-) | H S† |
| 10 | 0 | equal (+,+) | H |
| 10 | 1 | opposite (+,-) | H S† |
| 11 | 1 | equal (+,+) | H |
| 11 | 0 | opposite (+,-) | H S† |

Here S=diag(1,i). Let R be the composition of D and the conditional
re-encoding (d,a)->(d,t,r), and let C apply the six controlled block
transforms in the table, acting as identity on the four d=00 singleton states.
Then W=C R diagonalizes T. Indeed, R T R† exchanges r=0 and r=1 inside
each nonzero (d,t) block with the sign weights above, and C diagonalizes
each block. The d=00 states are fixed singletons with phases sigma_a. No final
R† is needed: the output computational ordering may be any eigenvalue order;
unrouting would only permute the output basis and incur extra cost.

## Alternative: common phase correction and d-controlled basis changes

The nonlinear conditional re-encoding above is not the only route. After the
two-CNOT difference transform R, each d sector acts on a as
T_d = CZ_A X_A^d, with X applied first and the CZ phase evaluated on the
flipped label. Apply the same diagonal correction
P_A = diag(1,1,1,-i) = CS† on A in every sector. Evaluating the phases on a
label a and its flipped label a xor d gives the corrected actions:

| d | Action of P_A T_d P_A† | Unconditional target basis change |
|---|---|---|
| 00 | CZ_A | none |
| 10 | phase i^(a1) times X0 | H on q0 |
| 01 | phase i^(a0) times X2 | H on q2 |
| 11 | phase i on even a0 xor a1, phase 1 on odd parity, times X0 X2 | H on q0 after CX(q0->q2) |

All displayed phase factors are diagonal in the computational labels. For the
d=11 block, CX(q0->q2) maps X0 X2 to X0 and maps the even/odd parity phase to
a diagonal phase on q2; H on q0 then diagonalizes X0. This gives an
alternative W = C P_A R, where C selects the indicated basis changes by d.

The controls can be arranged with one shared conditional parity change:
apply CX(q0->q2) only when d1=1; then apply H on q0 when d0=1 and H on q2
when d0=0,d1=1. The d=00 sector receives no H. In d=01, the conditional CX
commutes with the X2 flip and leaves its diagonal i^a0 phase intact; in d=11
it converts X0 X2 to X0 and the parity phase to a diagonal function of q2.
This is an exact by-hand branch rule, but its multi-control selectors still
need synthesis.

## Cost and next hypothesis

The difference transform D costs two CX gates abstractly. The exact Bell
decoder prefix can be realized as H_control followed by CX_control,target for
each pair, costing two CX total. Prefix plus D therefore consumes four of a
twelve-CX budget and leaves eight for the conditional router and basis
changes. In the shared-phase alternative, P_A costs two CX, leaving six for
its conditional CX/H selectors. These are arithmetic budgets, not
demonstrated decompositions. If the prefix's two native F gates are instead
counted at their generic two-CX compilation upper bound, prefix plus D costs
six CX and leaves six for the remaining transform (or four after P_A). The
full free-98 native family still has four CX plus four F gates and a generic
twelve-CX compilation upper bound. The fixed suffix has four native CX plus
two native F after its Bell prefix, and the q2-to-q3 difference route and
controlled selectors are not automatically present in that chronology.
Neither construction has been synthesized within the fixed suffix.

The actionable hypothesis is to synthesize this conditional router and
controlled block operation while omitting final unrouting, then compare its
actual native cost and chronology with the suffix. This is more targeted than
a fixed affine router and may guide an initializer if it fits. The affine
obstruction does not exclude arbitrary one-qubit rotations, native F
interactions, or the full 98-coordinate family. Run402's loss 0.375 and
invalid candidate remain one bounded numerical failure, not a family
exclusion or lower bound.
