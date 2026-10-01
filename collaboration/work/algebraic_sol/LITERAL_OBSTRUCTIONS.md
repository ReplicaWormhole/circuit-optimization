# Static audit of literal selector substitutions

Board 119; by-hand math and workflow artifact. No circuit matrix, symbolic
search, optimizer, or scientific-ledger computation was performed.
Qubit zero is most significant; gates are chronological and the target is
`U V4 = D U`. The accepted 13-CX circuit is unchanged.

Write `x=q0`, `y=q2`, `u=q1`, `v=q3` after the Bell decoder and difference
router. The routed shift is `T_d=CZ_xy X_x^u X_y^v`, with X acting first.

The literal chronological RCCX word

`H_t,T_t,CX(b,t),Tdg_t,CX(a,t),T_t,CX(b,t),Tdg_t,H_t`

has target blocks `I,I,Z,Y` for control values `00,01,10,11`.
All four blocks are Hermitian, so the complete RCCX is Hermitian and its
adjoint does not create a distinct variant. The original 16 choices contain
at most four different unitaries.

## CZ anchor obstruction

In sector `d=(u,v)=(1,0)`, the final wrapped CCH has both logical controls
zero, for either control order, and is identity. The first RCCX is identity
when its controls are `(v,x)`, or is `CZ_xy` when they are `(x,v)`.
With the substituted anchor `Btilde=CZ_xy`, the effective sector basis gate
is therefore `H_x CZ_xy` or `H_x`.

Restrict further to `y=1`. The original sector operator is
`T_10=Z_x X_x=iY_x`. Conjugating by Z and H only changes the sign of Y.
Thus the output still contains `+iY_x` or `-iY_x`, with off-diagonal
magnitude one. Every original literal choice fails exactly in this sector.
The independent verifier confirmed this by hand; see its board 121 report.

## Scope

This is a proof about these explicitly specified literal substitutions.
It excludes neither other relative-phase selectors, nor other 12-CX
circuits, nor the full98 family. No numerical screen is needed to establish
this narrow failure. No complete 12-CX proposal results from this audit.

Adversarial review checked chronological conjugation, the inactive branch,
and the distinction between adjoint and T-sign reversal. No blocker found
for the narrow obstruction; no gate-list simulation or exact checker ran.
