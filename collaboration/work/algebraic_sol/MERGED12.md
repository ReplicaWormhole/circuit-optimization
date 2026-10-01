# Four-CX phase-anchor/selector merge: complete structural 12-CX proposal

Board 123, continuing its bounded by-hand phase-repair derivation. No matrix
simulation, symbolic search, or optimizer was used. Independent full gate-list
numerical and exact checks remain required before any incumbent acceptance.
Keep q0 most significant, chronological gates, and `U V4 = D U`.

Use `x=q0`, `y=q2`, `u=q1`, `v=q3` after the specialized Bell decoder and
difference router. Start from the repaired 13-CX construction in
`PHASE_REPAIR.md`, but choose the first P RCCX controls `(x,v)`.
After conjugating its target by Sdg, P has branches:

- v=0: CZ_xy (inactive Z when x=1).
- v=1: CXxy (active X when x=1).

Therefore the chronological product of B=CSdg_xy then this P is a target-y
multiplexer with controls x,v and blocks:

| xv | U_xv = (P B)_xv |
|---|---|
| 00 | I |
| 01 | I |
| 10 | Z Sdg = S |
| 11 | X Sdg |

The P inactive phase is essential: choosing P controls `(v,x)` instead would
give a different block table and does not use this factorization.

## Four controlled reflections

Remove the control-only phase `T_x=diag(1,exp(i*pi/4))`. The determinant-one
target blocks become I, I, A, C, with

`A=exp(-i*pi/4) S = Rz(pi/2) = (I-iZ)/sqrt(2)`,

`C=exp(-i*pi/4) X Sdg = -i(X+Y)/sqrt(2)`.

Define target reflections

`N1=(-X+Y)/sqrt(2)`,

`N2=N4=(-X+Z)/sqrt(2)`,

`N3=-X`.

Chronologically apply controlled N1 with control x, controlled N2 with
control v, controlled N3 with control x, controlled N4 with control v,
then T_x. For each fixed control sector, the reflection product is

| xv | target product before T_x |
|---|---|
| 00 | I |
| 01 | N4 N2 = I |
| 10 | N3 N1 = (I-iZ)/sqrt(2) = A |
| 11 | N4 N3 N2 N1 = Z N1 = -i(X+Y)/sqrt(2) = C |

The last identity follows because reflection about
`(-X+Z)/sqrt(2)` maps `-X` to Z. All products are by-hand Pauli identities;
the output T_x restores the exact unstripped blocks, not just their
eigenvectors. This proves the four-CX merge at block level.

Every controlled reflection is locally equivalent to CX:

`N1=Rz(3*pi/4) X Rz(-3*pi/4)`,

`N2=N4=Ry(-3*pi/4) X Ry(3*pi/4)`.

Controlled N3 is CXxy accompanied by phase(pi) on x. Thus an explicit
chronological four-CX merged gate list is:

1. Rz_y(-3*pi/4)
2. CX(x,y)
3. Rz_y(3*pi/4)
4. Ry_y(3*pi/4)
5. CX(v,y)
6. Ry_y(-3*pi/4)
7. phase_x(pi)
8. CX(x,y)
9. Ry_y(3*pi/4)
10. CX(v,y)
11. Ry_y(-3*pi/4)
12. phase_x(pi/4)

Use `u3(theta=0,phi=0,lam=angle)` for both phase gates. The phase(pi) on
x implements the active minus sign in N3; dropping it changes the block
identity. The two control phase gates commute with every gate in this
merge, and may be combined as phase(5*pi/4) on x.

## Complete structural 12-CX circuit

Chronological blocks:

1. Bell decoder: `CX02,H0,CX13,H1` (2 CX).
2. Difference router: `CX01,CX23` (2 CX).
3. The explicit merged list above with x=0,y=2,v=3 (4 CX).
4. Exact CH(q1;q0):
   `Ry0(-pi/4),H0,CX10,H0,Ry0(pi/4)` (1 CX).
5. Final selector: negative control u=0, positive control v=1, target y.
   Apply `X1,Ry2(3*pi/4),Rx2(pi/2)`, then the literal RCCX with
   `(a,b,t)=(1,3,2)`, then `Rx2(-pi/2),Ry2(-3*pi/4),X1` (3 CX).

The final RCCX literal is chronological
`H2,T2,CX32,Tdg2,CX12,T2,CX32,Tdg2,H2`.
Its active Y becomes H under the stated target conjugation; its inactive Z
becomes Y and acts only on d=00. This inactive Y is a computational
permutation preserving diagonality in that sector.

Cost: `2+2+4+1+3=12 CX`, with no ancilla and rational-pi one-qubit gates.
No gates are assigned unknown or free entangling cost.

For the chosen P control order, after B/P and before the final rotations the
sector operators are CZ_xy, (-i)^y X_x, i^x X_y, and f_i(y) X_x for
d=00,10,01,11. CH turns the x exchanges into Z_x; final H_y turns the d=01
exchange into Z_y; final Y_y leaves the d=00 sector diagonal. This links the
exact merged block identity to the complete diagonalization derivation.

The resulting predicted diagonal is explicit:

- d=00: `(-1)^(x*(1+y))`.
- d=10: `(-i)^y (-1)^x`.
- d=01: `i^x (-1)^y`.
- d=11: `f_i(y) (-1)^x`, with `f_i(0)=i`, `f_i(1)=1`.

In computational order q0 q1 q2 q3 this is

`[1,1,1,-1,1,i,-i,1,-1,i,1,-i,-1,-i,i,-1]`.

Its root multiplicities are `(6,4,3,3)` for `(1,-1,i,-i)`. This predicted
order can be checked against the full candidate independently; arbitrary
output ordering remains allowed by the task.

## Review and limits

The independent verifier was asked to review the target axes, chronological
local conjugations, branch phases, and four block products before execution.
The coordinator owns reservation and all candidate checks. This note alone
does not change the incumbent or constitute independent gate-list acceptance.
Adversarial review checked the selector control order, phase(pi) sign,
T_x restoration, rational angles, and honest 12-CX total. No blocking issue
found in the by-hand derivation; serialized full-list and exact arithmetic
checks are still pending.
