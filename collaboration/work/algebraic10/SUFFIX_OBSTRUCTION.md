# Fixed-prefix joint suffix toward ten CNOTs

Board139, algebraic10. Artifact types: math and workflow. One bounded by-hand
investigation; no matrix evaluation, symbolic computation, optimizer, or
scientific search ran. The accepted exact11 circuit and its certificates
remain unchanged. Convention: q0 most significant, chronological gates,
`UV4=DU`. Write `x=q0,u=q1,y=q2,v=q3` after the Bell decoder/difference router.

## Exact joint suffix

The source is `experiments/runs/411/candidate.json`, with derivation
`collaboration/work/algebraic11nearhit/EXACT11_PROPOSAL.md`. Its final merged
reflection is controlled by v and acts on y:

`N'=-c X+s Z`, `s=sin(pi/8)`, `c=cos(pi/8)`.

It commutes through the following phase T on x and CH(u;x), since their
supports are disjoint. Thus, at the same gate cost, place N' immediately
before CZ(u;y) and the final controlled reflection

`B=s X+c Z`.

The complete three-CNOT joint target suffix has sector blocks, in (u,v):

| uv | F_uv |
|---|---|
| 00 | I |
| 10 | Z |
| 01 | B N' = -iY |
| 11 | B Z N' = A = (Z-X)/sqrt(2) |

These are exact Pauli identities. Using `XZ=-iY` and `ZX=iY`,

`B N'=-i(s^2+c^2)Y=-iY`,

`B Z N'=(s^2-c^2)X+2sc Z=(Z-X)/sqrt(2)`.

The first three blocks are computational monomial matrices (one nonzero
entry in each row and column). A has two nonzero entries in every row.
No sector scalar has been discarded: in particular F01 is -iY, not Y.

## What the prefix already diagonalizes

Let P be the original circuit through the first three merged reflections,
T_x, and CH(u;x), with N' moved into the suffix as above. The accepted
derivation gives the final sector operators

`D00=(-1)^(xy)`, `D10=(-i)^y Z_x`,
`D01=i^x Z_y`, `D11=f_i(y) Z_x`,

where `f_i(0)=i,f_i(1)=1`. Consequently its pre-suffix operators are
`Kuv=Fuv^dagger Duv Fuv`:

| uv | K_uv |
|---|---|
| 00 | CZ_xy |
| 10 | (-i)^y Z_x |
| 01 | -i^x Z_y |
| 11 | `[ (i+1) I/2 - (i-1) X_y/2 ] Z_x` |

Here `-i^x` means `-(i^x)`; it is not `(-i)^x`.
The identity `A Z A=-X` fixes the last sign. Thus, with this fixed prefix,
only uv=11 still needs a target basis change; the other three sectors already
have computational y eigenbases. In uv=00 the x=0 target is degenerate, but
the x=1 target is Z, so a target operation independent of x must preserve
the computational basis there. The uv=10 and uv=01 target eigenvalues are
distinct. All statements in this section are derived from the certified
source sector table, not a new full gate-list certificate.

## Narrow two-interaction obstruction

Freeze P. Restrict the replacement suffix to sector-preserving, target-y
controlled gates with u and v each used once. Allow arbitrary unconditional
target rotations before, between, and after those gates, and arbitrary
one-qubit diagonal phases on u and v. For chronology u then v, its blocks
can be written

`Guv=L2 Nv^v L1 Nu^u L0`.

Then, with no commutation assumption on the target rotations,

`G11=G01 G00^-1 G10`.

The reverse chronology gives `G11=G10 G00^-1 G01`. One-qubit control phase
factors cancel in these relations. The relation holds even when Nu and Nv
are arbitrary target unitaries rather than reflections; hence the
obstruction does not depend on a special choice of CX conjugating axes.

Allow sectorwise output y flips and phases, represented by arbitrary target
monomial M_uv. A replacement with the same fixed x,u,v labels must have
`Guv=Muv Fuv`: target eigenvalues in uv=10 and uv=01 are nondegenerate, and
uv=00 at x=1 is nondegenerate. A y operation independent of x must therefore
be monomial in these three sectors. The desired uv=11 eigenspaces require
`G11=M11 A`; multiplying A on the left by a monomial matrix cannot remove
either nonzero entry in any row. But the two-interaction relation makes
G11 a product of three monomial matrices, hence monomial. Contradiction.

Therefore there is no saving from three to two target interactions within
this fixed-prefix, sector-preserving template, including its conditional
target basis freedom and its permitted output target permutations/phases.
The two-CX product would require its fourth block to be monomial, while
diagonalizing the remaining X target axis requires full target support.

## Remaining scope and actionable next hypothesis

This does not exclude a ten-CNOT circuit. It does not cover changing any of
the first three merged reflections, the Bell decoder/router, or CH(u;x),
nor suffixes that mix or permute x,u,v labels, nor changes to a degenerate
eigenbasis that mix different sectors. In particular, an arbitrary global
output computational permutation is broader than the sectorwise y
permutations analyzed here. Two-CX suffix circuits with changed control
bases are also outside this argument.

The next distinct hypothesis is to free the earlier merged reflections
jointly with CH and the final selectors; a local deletion refit may find
new sector axes that avoid the single exceptional uv=11 corner. The
numerical collaborator's fresh exact11 deletion preparation already tests
this broader possibility with fully free one-qubit layers. A near-hit would
identify which prefix interaction/basis must change and justify a new
by-hand construction. The present note is not a topology exclusion for
that fit and gives no resource authorization to execute it.

## Independent exchange and adversarial review

The reduced suffix blocks and product-law obstruction were sent directly to
numerical10 and the independent verifier (live agent numerical_sol, board
owner verifier10), with persisted board139 handoffs. The verifier
independently confirmed both Pauli products, the product law, and its
monomial contradiction by hand. No new ten-CNOT word or predicted diagonal
is proposed. The incumbent eleven-CNOT diagonal is preserved.

Review followed `~/.codex/adversarial_review.md`. Checked chronological
order, commuting supports, XZ/ZX signs, the sector scalar -i, the x=0
degeneracy in uv=00, the target support argument, control-phase cancellation,
and the restricted output-permutation scope. No blocking issue was found
within the stated template. The broader ten-CNOT problem remains open and
no matrix or proof-assistant validation was performed for this note.
