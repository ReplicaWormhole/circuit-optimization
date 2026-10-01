# Independent output-routing freedom audit

Board130. This is by-hand derivation and source-note review only; no matrices,
symbolic search, objective, optimizer or exact checker has run. The certified12
candidate is unchanged. The argument below restricts literal/router-preserving
changes and does not exclude the forthcoming fully-local deletion refits.

Let E be the Bell decoder, R=R01 R23 the difference router, and F the entire
following tail (mergedBP, CH, finalselector). The accepted word is U=F R E.
Qubit0 is most significant; x=q0,y=q2,u=q1,v=q3 are the post-router labels.
R01=CX(x->u), R23=CX(y->v). Every completed tail block preserves u,v: mergedBP
acts on y with controls x,v; CH acts on x with control u; the final selector's
X_u wrappers cancel and it acts solely on y conditioned on u,v.

## Why literal omission cannot be repaired by output permutations

Omit R01 and retain R23. In Bell coordinates, the original shift acts
(a,b)->(b,a), with its known scalar Bell sign. After R23 alone, u is the
old b0 while x is old a0. The transformed shift maps u to x; for x!=u it
therefore has a nonzero block between distinct u sectors. F is block diagonal
in u, so F conjugation cannot make that block zero: unitary multiplication
within the two blocks preserves whether their connecting operator vanishes.
A computational diagonal operator has no such block. Thus literal omission
of R01 fails, as does any replacement tail that still preserves u.

Omit R23 and retain R01. Analogously v is old b1 while y is old a1, and the
transformed shift maps v to y. It has nonzero connecting blocks for y!=v.
The tail preserving v cannot eliminate these blocks. Literal omission of
R23 fails, as does any replacement tail that still preserves v.

An output computational permutation or diagonal phase reindexes/phases
existing off-diagonal entries; it cannot turn a non-diagonal matrix into a
diagonal one. Merely renaming the output sector labels is likewise not a
physical replacement for the missing coherent router operation. Arbitrary
orthonormal changes within degenerate eigenspaces do not invalidate the block
argument: if the resulting target were diagonal its intersector blocks would
have to vanish. A synthesis that absorbs a router while changing the tail's
control structure remains possible and is not excluded by this argument.

## Independent explicit output-residue checks

These confirm algebraic_sol's two by-hand examples. Write M for mergedBP,
C=CH and Hfin for the final relative-phase selector. Then F=Hfin C M.
Removing R01 changes the output by
Q01=Hfin C R01 Cdag Hfin-dag, since R01 commutes with M.
On input |x,u,y,v>=|0,1,0,0>, Hfin-dag is I; Cdag produces the equal x
superposition. R01 sends its x=1 branch to u=0. C produces two amplitudes
1/2 in the u=1 branch, while Hfin=Y_y in the u=0,v=0 branch gives

Q01|0100> = (1/2)|0100> + (1/2)|1100> + (i/sqrt2)|1010>.

Thus Q01 is not an output computational permutation up to basis phases.
For R23 the residue is
Q23=Hfin C M R23 Mdag Cdag Hfin-dag. On |0000>, Hfin-dag=Y_y gives
i|0010>; C and Mdag are I along that branch. R23 flips v to1, M and C
remain I, and Hfin=H_y gives

Q23|0000> = (i/sqrt2)(|0001>-|0011>).

Again the residue is not monomial. These examples independently confirm
narrow literal-deletion failure; a nonmonomial output residue alone would
not prove failure if it mixed only equal eigenvalues, which is why the
intersector target-block argument above is the actual diagonalization audit.

## Consequence for the bounded fit

All144 single-qubit rotation-vector coordinates in an11CX deletion word may
change the physical control basis and break the tail's u/v preservation.
The proposed12-case fully-local fit is therefore a distinct bounded numerical
hypothesis; the fixed-tail obstruction does not exclude it. A promising fit
needs complete gate-list independent checks and a new exact certificate
before changing the accepted12 result. Failed fits provide no global bound.

Adversarial review used ~/.codex/adversarial_review.md. It checked both missing
router conjugations, phase/control order, residue amplitudes, fixed-tail versus
fully-free scope and output-permutation claims. No unresolved issue found.
This directory is intentional research evidence; concurrent user work preserved.
