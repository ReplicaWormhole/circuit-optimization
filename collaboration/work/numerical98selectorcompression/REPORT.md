# Phase-aware compression of Bell-difference selectors

Board 115. Read-only synthesis analysis based on the common-phase route in
algebraic98difference. No matrices, gate decompositions, circuit simulations,
optimizer, or search were executed here.

## Fixed target and cost arithmetic

The exact abstract diagonalizer after the Bell prefix is C P_A R:

- R: CX(q0->q1) and CX(q2->q3), two CX.
- P_A: CS† on (q0,q2), two CX in an exact controlled-phase decomposition.
- C: a Toffoli selector (controls q3,q0; target q2), a controlled-H selector
  (q1 controls H on q0), and a doubly controlled-H selector (q1=0,q3=1
  controls H on q2).

The specialized exact Bell prefix costs two CX, so prefix + R + P_A costs
six CX and leaves six in a 12-CX budget for C. The literal standard selector
implementations cost 6 + 1 + 6 = 13 CX, for total 19. Replacing both 6-CX
selectors by 3-CX relative-phase candidates would reduce the arithmetic to
13, but does not itself preserve the diagonalizer. One more CX must be saved
by a phase-aware merge with P_A or another selector, and the proposed CCH
candidate must be an actual decomposition before this count has meaning.

The generic native-family count remains four CX plus four F gates, with a
12-CX compilation upper bound. The two native F gates in the Bell prefix can
be specialized to the exact two-CX Bell decoder, but that does not change the
unsearched family's topology or prove the remaining selectors fit.

## When relative phases can be harmless

Let U=C P_A R be the exact basis change and D=U T U† be diagonal. Replacing U
by Q U is safe if Q is diagonal in the final computational output basis:
Q D Q†=D. More generally, a diagonal phase inserted in an intermediate
Bell-label basis is harmless only if it is constant on each exchange block
(equivalently, commutes with the routed operator there), or if a later exact
correction moves its entire effect to the final output side. A phase that
depends on the orientation within a pair is not block-constant.

This distinction matters for the two candidate savings:

1. A 3-CX relative-phase Toffoli may be written as the exact Toffoli times
   diagonal phase factors, but their placement and dependence on q0,q2,q3
   must be retained. The Toffoli occurs before the controlled Hadamards.
   Its residual phase is not automatically a final output phase and can fail
   to commute through a later H on q0 or q2.
2. A doubly controlled H is not a Toffoli. If a candidate is formed by
   target-conjugating a relative-phase Toffoli, its residual correction is
   conjugated by that target basis change: V Delta V†. A diagonal Delta on
   the target generally becomes non-diagonal, so it cannot be discarded as a
   diagonal output phase. In inactive-control branches, q2 still labels
   eigenvalue data (for example d=00 has a CZ-dependent phase, d=10 has an
   i^(a1) phase, and d=11 has an i^(1-a1) phase after the shared selector).
   Thus target-dependent relative phases can mix distinct roots even where
   the intended controlled H is inactive.

The shared P_A=CS† phase is a plausible place to absorb a residual phase only
after the candidate decompositions expose their exact diagonal phase
polynomials and those polynomials are checked against the branch eigenvalue
labels. No such cancellation is established here, and P_A cannot simply be
omitted on the current derivation.

## Recommended finite next step

Before a new fit, freeze literal exact and relative-phase gate lists for the
two selectors, including all local phase corrections and control polarities.
Use a finite budget of at most 16 variants: the exact baseline, then
predeclared forward/inverse relative-phase choices for each selector and
predeclared CS†/CS sharing choices, with any remaining slots assigned only
to explicit phase-correction placements. For each variant, compare the full
four-qubit unitary against C P_A R up to global phase and verify that the
resulting conjugated target is diagonal with the required root multiplicities.
Only variants that pass that complete-unitary check may be costed as
diagonalizers. The test must include inactive-control branches; a truth table
for the classical Toffoli action alone is insufficient.

If no literal 12-CX variant survives, the alternative is one separately
reserved bounded all-free 12-CX fit with a newly frozen initialization and
budget. That would be numerical evidence only. Neither the 19-CX construction
nor failure of this finite variant screen proves a 12-CX lower bound or
excludes the full family.
