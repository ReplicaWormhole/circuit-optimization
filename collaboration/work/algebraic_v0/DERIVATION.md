# Phase-compatible reduced v=0 old and reordered sources

Board155, algebraic_v0. Artifact types: math, computation specification and
workflow. Static by-hand derivation only: no matrix evaluation, symbolic
search, RNG, fitted-angle interpretation or optimizer executed here.
Chronological gates; q0 is MSB; n3 wires x=q0,u=q1,y=q2. The reduced target
is `W=CZ_xy CX(u,x)`, with CX acting first. This is not the cyclic shift
on three physical qubits. The source files require a checker of this
explicit reduced target, not the default n3 cyclic-shift target.

## Independent Bell-label phase derivation

Use opposite-pair Bell labels a,b on physical wires0,2 and c,d on1,3,
with pair state `(|0,b>+(-1)^a |1,1 xor b>)/sqrt(2)` and its c,d analogue.
Its physical computational terms have bits

`(j,l,j xor b,l xor d)` and coefficient `(-1)^(a*j+c*l)/2`.

The right shift sends these bits to
`(l xor d,j,l,j xor b)`. Writing the new first pair bit as
`j'=l xor d` and new second pair bit as `l'=j`, its phase becomes
`(-1)^(c*d) (-1)^(c*j'+a*l')`. Thus the exact decoded shift acts as

`|a,c,b,d> -> (-1)^(c*d) |c,a,d,b>`.

The accepted decoder is chronological CX02,H0,CX13,H1. Its difference
router CX01,CX23 makes `x=a,u=a xor c,y=b,v=b xor d`. The output labels are
`x'=x xor u,u'=u,y'=y xor v,v'=v`; the phase is
`(-1)^((x xor u)*(y xor v))=(-1)^(x'*y')`. Hence the routed operator is
exactly `CZ_xy X_x^u X_y^v`, including every phase. At v=0:

`W00=CZ_xy`, `W10=[[0,I_y],[Z_y,0]]` in x blocks.

At fixed (u,y), these x operators are respectively
`I,Z,X,iY` for 00,01,10,11. In particular the last is iY, not Y or -iY.

## Literal old source

`old_source.json` is the accepted411 v=0 restriction with inactive v
reflection wrappers removed in exact inverse pairs. Its chronology is

1. `Rz_y(-3pi/4),CXxy,Rz_y(3pi/4)`.
2. `Rz_y(-pi),CXxy,Rz_y(pi)`.
3. `U3_x(0,0,pi/4)`.
4. `Ry_x(-pi/4),H_x,CXux,H_x,Ry_x(pi/4)`.
5. `H_y,CXuy,H_y`.

The first target reflections are `N1=(-X+Y)/sqrt(2)` and `N3=-X`.
Their active product is `N3 N1=(I-iZ)/sqrt(2)=Rz(pi/2)`; multiplying
the active phase exp(i*pi/4) on x gives S_y=diag(1,i), with inactive
x=0 block I. Therefore its first three blocks realize controlled-x S_y
exactly, including the active minus sign and the full U3 phase.
The fourth block is exact CH(u;x). The last CZ(u;y) commutes with the
resulting diagonal operator. The old predicted D in MSB x,u,y order is

`[1,1,1,-i,1,-1,-1,i]`.

This source has15 elementary gates and4 CX, topology xy,xy,ux,uy, and no
one-qubit u gates. Controlled S followed by CH already diagonalizes the
reduced target at3 CX; the final CZ is a commuting output phase included
to match the requested four-interaction extension template. No claim that
four CX is the minimum reduced cost is made.

## Reordered source with tilted x basis

Define exact x reflections

`K=(X+Y)/sqrt(2)`, `J=(Y+Z)/sqrt(2)`,
`D=(K+Z)/sqrt(2)=(X+Y)/2+Z/sqrt(2)`,
`K0=Z D Z`, `R=Z J Z=(Z-Y)/sqrt(2)`.

K and Z are orthogonal Pauli axes, so D is a unit reflection and
`D Z D=K`. Chronologically apply

`D_x, CZ_xy, K0_x, controlled-R(u;x), CZ_xy, CZ_uy`.

This has the requested CX order xy,ux,xy,uy. Computational u is preserved;
there are no one-qubit u operations. Before final CZ_uy, its x blocks in
fixed (u,y) are

| uy | U_uy |
|---|---|
| 00 | K0 D = ZK |
| 01 | Z K0 Z D = I |
| 10 | R K0 D = ZJK |
| 11 | Z R K0 Z D = J |

The final CZ_uy multiplies the last block by -1 and leaves all others
unchanged. This sign is retained in the complete word. It does not change
operator conjugation. Since `K X K=Y` and `J Y J=Z`, the four target x
operators I,Z,X,iY become I,Z,Z,iZ. Thus the new predicted D is

`[1,1,1,i,1,-1,-1,-i]`.

The roots i and -i exchange output positions relative to the old source;
this is permitted output eigenvalue ordering, not a lost phase. Both lists
have multiplicities4,2,1,1 for1,-1,i,-i.

Every reflection has rational-pi one-qubit wrappers:

`D=Rz(pi/4) Ry(pi/4) Z Ry(-pi/4) Rz(-pi/4)`,

`K0=Rz(-3pi/4) Ry(pi/4) Z Ry(-pi/4) Rz(3pi/4)`,

`R=Rx(pi/4) Z Rx(-pi/4)`.

The exact Z in the first two lines is U3(0,0,pi); replacing it by Rz(pi)
would change the phase. Controlled R is chronological
`Rx_x(-pi/4),H_x,CXux,H_x,Rx_x(pi/4)` and has inactive block I exactly.
`new_source.json` spells out all24 gates. Its two adjacent final H_y gates
are deliberately retained, so no wrapper cancellation must be inferred.
The initial D and subsequent K0 rotate x non-diagonally, violating the
computational-x assumptions of the earlier early-CH obstruction.

## Provenance and remaining full problem

SHA256 hashes, computed only as file bookkeeping:

- old_source.json: `93030d184db2725c2b32e68b8b28eb8481727abbdf9f38af9018771b32842b82`.
- new_source.json: `46a3b51463e88207caf735114d8f8331b8a1ed2b2eec34debab4c015e565cd3c`.

The independent verifier confirmed the Bell phase, old-source N3N1/T
repair, all new sector blocks, the D/R wrapper signs and predicted roots
by hand. Complete saved-list numerical and exact certification must still
precede a computationally certified reduced result. Root owns reservation
and execution; this finding can supersede the prepared reduced fit rather
than execute an unnecessary optimizer.

The proposed full ten-CX comparator prepends the accepted exact four-CX
decoder/router and inserts a bare v-to-y CX after the first xy interaction
and another at the end. For this literal new source, independent by-hand
review identifies a definite v=1 failure: the first insertion lies inside
the H_y,CXxy,H_y wrapper and hence applies Z_y in v=1; the final insertion
applies X_y. All net reduced blocks commute with Z_y, so the full v=1
suffix is exactly `X_y Z_y Ured`. The routed operator
`Wv1=CZ_xy X_x^u X_y` anticommutes with Z_y. Ured's conjugation preserves
that anticommutation, and the final monomial X_y Z_y still preserves it
(it sends Z_y to -Z_y). Thus the resulting v=1 operator cannot be
computationally diagonal. This applies only to the specified literal
bare-insertion comparator, not the fully free topology.

No full V4 upper bound, acceptance, or incumbent change follows. The
precise next constructive hypothesis is to replace the two bare v
insertions by adjustable controlled target reflections and derive both
v=1 sectors jointly while preserving the now-solved v=0 blocks. It must
include a non-computational target-y basis in v=1, rather than only
monomial target operations. This has no separate computation authorization.

For that distinct general-reflection extension, the first controlled
reflection is inserted after the complete first CZxy wrapper: after
primitive8, zero-based source position7 (H_y), and before source
position8 (Rz_x(3*pi/4)). The final reflection is after all24 primitives,
after zero-based position23 (H_y). These cuts apply to the unchanged
new_source.json hash above and yield suffix CX order02,32,10,02,12,32.
The first cut is outside the wrapper, unlike the known-failed literal
bare comparator. General reflection axes at these new cuts are untested.

## Adversarial review

Review followed ~/.codex/adversarial_review.md. Checked Bell encoder
conventions, right-shift order, router input/output phase, W10's lower Z,
reflection products, chronological rotation signs, exact Z versus Rz,
inactive controlled-R identity, final sector scalar -1, root order,
both gate counts, u preservation, the literal bare-v anticommutation
obstruction, and reduced-versus-full scope. No
blocking issue found by hand. Scripts/checkers, matrices, or searches
were not executed by this collaborator. New source JSON and this note
are intentional research evidence; all prior artifacts remain untouched.
