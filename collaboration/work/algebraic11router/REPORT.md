# Difference-router absorption after the certified 12-CX circuit

Board 129, algebraic11router. One bounded by-hand investigation only: no
matrix evaluation, symbolic search, optimizer, or new scientific run.
Read accepted run406 and `collaboration/EXACT12_ACCEPTANCE.json`; the
certified12 candidate and all previous evidence remain unchanged.

Use x=q0,u=q1,y=q2,v=q3 after the Bell decoder. The chronological12 word is
E, R01, R23, M, CH, F. Here R01=CX(x,u), R23=CX(y,v), M is the exact4CX
phase-anchor/selector merge, CH is controlled-H(u;x), and F is the repaired
final selector. In routed coordinates, F has branches Y_y for d00, H_y for
d01, and I for d10/d11. The target remains UV4=DU, q0 most significant,
arbitrary output root order and arbitrary degenerate eigenbases.

## A useful exact rearrangement

R01 commutes with R23 and with M: M acts on y, preserves x,v, and ignores u.
Therefore R01 can move immediately before CH without changing the unitary:

`E, R23, M, R01, CH, F`.

This exposes a two-CX x/u synthesis target `G=CH(u;x) CX(x,u)`. It is a
changed architecture placement, rather than a smaller circuit.

On input basis states of (x,u), G has the exact columns

- 00 -> 00.
- 01 -> (01+11)/sqrt(2).
- 10 -> (01-11)/sqrt(2).
- 11 -> 10.

No one-CX realization up to a permitted full-circuit eigenbasis change was
derived. Replacing CH by an unconditional H saves one CX but breaks the
d00/d01 sectors, which have genuine x-dependent roots. Degeneracy in some
sectors does not authorize this basis change in all sectors.

## Why literal deletion and a sector-preserving tail fail

The existing M, CH, F tail preserves u and v exactly. The X_u pair wrapping
F restores u; its complete operation is controlled by u,v and changes y.

If R01 is omitted but R23 retained, after the Bell-coordinate exchange the
unrouted u bit maps to x. For x!=u, the partially routed shift has a nonzero
block between distinct u sectors. An arbitrary invertible tail preserving u
conjugates such a block as `C_u A C_uprime_dagger`; it cannot turn a nonzero
A into zero. Therefore no u-preserving modified reflection/selector tail can
diagonalize the partially routed shift. This includes arbitrary bases inside
its u sectors, not just the literal original gates.

Similarly, omitting R23 with R01 retained leaves v mapping to y. States y!=v
give a nonzero off-v block, which no v-preserving tail can remove. An output
computational permutation reindexes both matrix indices and cannot turn a
nonzero off-diagonal entry into a diagonal one.

Concrete checks on the literal source tail give nonmonomial deletion errors.
Deleting R01 gives, on the output basis state |x0,u1,y0,v0>,

`(1/2)|0,1,0,0>+(1/2)|1,1,0,0>+(i/sqrt(2))|1,0,1,0>`.

Deleting R23 gives on |0,0,0,0>,

`(i/sqrt(2))(|0,0,0,1>-|0,0,1,1>)`.

These checks rule out repairing the exact omitted-router circuit solely by
an output permutation. They do not rule out a joint synthesis that changes
the spectator label, or fully free local-gate refits after a deletion.
The independent verifier confirmed the sector-preservation obstruction.

## Actual joint R23 rewrite

The x-controlled blocks of M are particularly simple. From the certified
first-selector order (x,v),

`M_x0=I`,

`M_x1=CX(v,y) S_y CZ_yv`.

This follows from target blocks S_y when v=0 and X_y Sdg_y when v=1.
Put R=R23. Then the exact joint conjugation R M R has x=0 block I and x=1
block

`J=SWAP_yv CZ_yv Sdg_y`.

Indeed `R CX(v,y) R=SWAP_yv`, R commutes with S_y, and
`R CZ_yv R=Z_y CZ_yv`; combine S_y Z_y=Sdg_y. No SWAP is being assigned
zero physical cost: this is an exact algebraic synthesis target.

J acts on y/v basis as

- 00 -> 00.
- 01 -> 10.
- 10 -> -i 01.
- 11 -> i 11.

Its even block diag(1,i) and odd block [[0,-i],[1,0]] both have determinant
i. Thus J is a parity-preserving, phase-weighted fermionic swap, in the
explicit determinant-matched two-qubit sense above. This observation is a
structural identity, not a claim about its controlled CNOT synthesis cost.

Allowing output relabeling by R gives a complete transformed synthesis
problem. Since CH commutes with R,

`R U = (R F R) CH (R M R) R01 E`.

The original R23+M joint cost is5CX. An11CX construction would result if
the controlled-x J could be built with4CX AND the conjugated final selector
Fprime=R F R still cost3CX. Neither cost is proved here. The transformed
final selector changes v, so it evades the earlier sector-preservation
obstruction; it cannot simply be left equal to F.

More explicitly, for u=1 Fprime=I. For u=0 it acts according to parity
p=y xor v: it is the conjugated Y block at p=0 and conjugated H block at
p=1. Using Pauli conjugation by CX(y,v),

`Fprime_u0 = [(I+Z_y Z_v)/2] Y_y X_v`
`             + [(I-Z_y Z_v)/2](X_y X_v+Z_y)/sqrt(2)`.

The displayed projectors commute with their corresponding target operators.
This is an explicit joint target for a finite future synthesis experiment;
it is not a gate-list candidate with known11CX cost.

The independent verifier confirmed the joint J identity and gave the simpler
parity-block form of Fprime_u0: even states (00,11) carry Y, while odd states
(01,10) carry H. Both blocks have determinant -1. Thus both J and Fprime_u0
are determinant-matched parity-preserving two-qubit operators. The control-x
and control-(u=0) synthesis costs remain unproved; the matchgate designation
does not assign a CNOT cost.

The verifier/root also identified a narrow failure of a proposed3CX template
for Fprime: conjugating a single controlled local Pauli by a one-CX pair basis
change would require operator Schmidt rank at most two for Fprime_u0.
Its Pauli coefficient form is
`(XY+YX)/2 + (XX+YY+ZI-IZ)/(2*sqrt(2))`, which has rank four: the X/Y
coefficient block has determinant -1/8 and ZI/IZ contribute two independent
entries. This rules out that simple template, not general3CX controlled
matchgate synthesis. This externally supplied check is recorded as a peer
finding; it was not independently re-derived in this bounded investigation.

## Actionable scope and closeout

- For R01 absorption, a changed joint merge must allow u to change. Moving
  R01 next to CH exposes the relevant x/u two-gate block.
- For R23 absorption, a changed joint merge must allow v to change. The exact
  controlled-J target and Fprime above express that freedom coherently.
- Fully free local deletion refits may leave the old preserved-register class
  and are not excluded by this report.

No explicit11CX circuit was found within this bounded derivation. No
optimality, family exclusion, or new lower bound follows. Adversarial review
checked chronological matrix order, omitted-router block reasoning, phases
in M_x1, SWAP algebra, and the need to conjugate the final selector. No
blocking issue found in these narrow identities. No numerical/exact matrix
check ran; numerical cost/implementation claims remain open.
