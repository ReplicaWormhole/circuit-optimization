# Early CH needs a different x basis

Board149, algebraic10reorder. Bounded by-hand structural work only; no
executed matrix evaluation, randomness, symbolic search, angle recognition,
or optimizer. Preserve all exact11/12/13/14 certificates. Qubit0 is MSB;
gates are chronological; the target is UV4=DU. Set x=q0,u=q1,y=q2,v=q3.

The proposed complete CX word is

`02,13,01,23 | 02,32,10,02,12,32`.

Its exact decoder/router prefix costs four CX and its reordered suffix six,
giving ten total. Moving the u-to-x interaction before the second x-to-y
interaction changes the word. The bare chronological triangle
`CX02,CX10,CX12 = CX10,CX02` is correct, but the incumbent's one-qubit
reflection/CH wrappers prevent applying it as a cancellation identity.
No complete rational-pi ten-CX candidate was found in this investigation.

## The quarter-turn phase before early CH

Undo the known certified12 merged block and CH in uv=10. The merged block
is `M=diag_x(I_y,S_y)`, where `S=diag(1,i)`; after CH the sector is
`D10=Sdg_y Z_x`. Consequently the exact decoder/router sector is

`W10=Mdag (Sdg_y X_x) M = [[0,I_y],[Z_y,0]]`.

These are x blocks: the upper off-diagonal target block is `Sdg S=I`,
the lower is `Sdg^2=Z`. Keeping this lower Z retains the relative
quarter-turn phase; replacing W10 by a simple X_x exchange is invalid.
Equivalently, target y=0 has x operator X and y=1 has x operator iY.

## A restricted obstruction for moving CH earlier

Assume u,v remain computational sector labels. Before early CH, x is
still a computational control, and the sole x-to-y interaction has the
form of one controlled reflection, allowing arbitrary target-local
wrappers and diagonal x phases. Thus its sector blocks are Q0,Q1 with

`Q=diag_x(Q0,Q1)`, `E=Q1 Q0dag=exp(i phi) N`, `N^2=I`.

In particular `E^2` is scalar. The v-to-y interaction is inactive in
v=0; its unconditional target wrappers are included in Q0,Q1. Conjugation
of W10 gives

`Q W10 Qdag = [[0,Edag],[E P,0]]`, `P=Q0 Z Q0dag`.

P is a non-scalar target reflection with eigenvalues +1,-1. Let the early
u-to-x gate act by any x-only unitary K in u=1; this allows more freedom
than the fixed exact CH wrapper. For K to make the displayed operator
computational-x block diagonal, it must commute before that rotation with
`Kdag Z_x K`. Write that axis as `k_x X+k_y Y+k_z Z`. Its commutator's
off-diagonal x blocks force k_z=0, since E is invertible. The diagonal
x blocks then require the two target blocks E P and Edag to be
proportional: `E P=rho Edag`, with scalar rho. But this gives

`P=rho (Edag)^2`,

which is scalar because E^2 is scalar. Contradiction. Therefore no x-only
early CH can diagonalize this sector after just one computational-x
controlled reflection and target-local wrappers.

If every subsequent gate preserves computational x, its action is
block diagonal in x. It sends a nonzero off-diagonal target block A to
`L0 A L1dag`, which cannot vanish for unitary L0,L1. In particular the
remaining computational-x controlled reflection and the u/v controlled
target gates cannot repair this failure without further x basis changes.
Diagonal x phases, sector scalar phases, and output target flips do not
remove the contradiction.

## Scope and constructive escape

This excludes transplanting an early CH into that restricted reflection
template. It does not exclude the fully free 132-coordinate reordered
family. A non-diagonal x rotation before the first controlled reflection
invalidates the Q assumption; a later non-diagonal x rotation invalidates
the subsequent block-diagonal restriction. Mixing u/v sector bases also
escapes it. No lower bound, topology exclusion, or exact interpretation of
run416's failed angles follows.

The next constructive algebraic hypothesis is to solve the two v=0 sectors
jointly with tilted x control axes and a final x readout rotation after
the second x-to-y interaction. Require their complete operators W00=CZ_xy
and W10 above to become diagonal, retaining the y-dependent quarter-turn
phase. This asks the two x-to-y interactions to synthesize the phase anchor
across the early u-to-x gate, rather than freezing their control axes and
moving the old CH wrapper literally. After obtaining a phase-compatible
v=0 solution, test whether the two v-to-y reflections can extend it to both
v=1 sectors. This is a by-hand subproblem already within the planned full
local topology's freedoms, not a request to repeat an unchanged deletion
or authorize another numerical run. No axes or angles are yet asserted to
solve it.

## Review and handoffs

The independent verifier (numerical_sol, owner verifier10reorder)
confirmed W10's phases, Q's blocks, the proportionality condition, and the
scalar-square contradiction by hand. Findings were sent directly to
numerical10 and the coordinator and persisted as board149 handoffs.
Adversarial review followed ~/.codex/adversarial_review.md, checking the
chronology, exact resource count, x-block convention, Sdg signs, invertible
blocks, the later-x restriction, and claims versus evidence. No blocking
issue was found within the stated assumptions. General reordered-family
existence remains open; full gate-list and exact arithmetic checks would
be required for any actual candidate.
