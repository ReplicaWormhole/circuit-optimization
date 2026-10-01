# Joint merged-block/final-selector eleven-CX proposal

Board 135, algebraic11nearhit. This is a by-hand structural derivation from
the certified12 word and run409's deletion9 topology. No matrix evaluation,
symbolic/numerical angle-recognition search, or optimizer ran here. The
proposal needs independently reviewed complete gate-list numerical and exact
certification; it does not yet change the incumbent.

Convention: chronological gates, q0 most significant, UV4=DU, no ancilla.
Use x=q0,u=q1,y=q2,v=q3 after the Bell decoder/difference router. The exact12
source is `experiments/runs/406/candidate.json`. Its merged block M has four
controlled target-y reflections, chronological N1_x,N2_v,N3_x,N4_v, followed
by T_x:

`N1=(-X+Y)/sqrt(2)`, `N2=N4=(Z-X)/sqrt(2)`, `N3=-X`.

## Rotate the merged output basis without another CX

Let L=Ry_y(-pi/4), so `L X Ldag=H=(X+Z)/sqrt(2)` and
`L Z Ldag=A=(Z-X)/sqrt(2)`. Replace only the fourth reflection by

`N4prime=L N4=-cos(pi/8) X+sin(pi/8) Z`.

This is Hermitian and squares to I: L rotates about Y, perpendicular to the
N4 axis in the X/Z plane, and the exact product below establishes the
Hermitian reflection. Explicitly,

`N4prime=Ry(-7*pi/8) X Ry(7*pi/8)`.

The new fourth controlled reflection still costs one CX. In v=0 nothing
changes; in v=1 its product is L times the old fourth reflection. Therefore
the exact new merged block is `Mprime=C_v(L_y) M` at the same four-CX cost.
This uses the existing final reflection to absorb a conditional quarter-turn;
it does not separately compile C_v(L_y).

## Two-CX final selector

After old M and exact CH(u;x), the sector operators were

| uv | operator |
|---|---|
| 00 | CZ_xy |
| 10 | (-i)^y Z_x |
| 01 | i^x X_y |
| 11 | f_i(y) Z_x, f_i(0)=i, f_i(1)=1 |

The new conditional L rotates X_y to H_y and Z_y to A_y on v=1. Thus the
u0v1 branch needs H_y diagonalized, while the u1v1 branch needs A_y
diagonalized. All v=0 branches are already diagonal on y.

Choose `B=sin(pi/8) X+cos(pi/8) Z`. Reflection B maps H to Z, while
reflection Z maps A to H:

`B H B=Z`, `Z A Z=H`.

Apply chronologically controlled-Z(u;y), followed by controlled-B(v;y).
These are two one-CX gates. Their branch actions are

| uv | final action | resulting target basis |
|---|---|---|
| 00 | I | old diagonal CZ unchanged |
| 10 | Z_y | old diagonal phase unchanged |
| 01 | B_y | H_y becomes Z_y |
| 11 | B_y Z_y | A_y first becomes H_y, then Z_y |

All four sectors therefore become diagonal. In particular, a two-CX final
selector is possible after changing the merged output basis. It does not
realize the former fixed final selector by itself and contradicts none of
the earlier narrow fixed-template obstructions. No control-register basis
change is needed for this explicit joint construction.

The d11 phase is preserved as an operator, including its scalar part:
`f_i(y)=(i+1)I/2+(i-1)Z_y/2`. The conditional L maps only Z to A; the
final Z and B map A through H back to Z. Thus there is no dropped sector
global phase or implicit change of the i/-i roots.

## Complete gate chronology

1. Bell decoder `CX02,H0,CX13,H1` (2 CX).
2. Difference router `CX01,CX23` (2 CX).
3. Modified merged block (4 CX), chronological:
   `Rz2(-3*pi/4),CX02,Rz2(3*pi/4),`
   `Ry2(3*pi/4),CX32,Ry2(-3*pi/4),`
   `Rz2(-pi),CX02,Rz2(pi),`
   `Ry2(7*pi/8),CX32,Ry2(-7*pi/8),T0`.
4. Exact CH(u;x), chronological
   `Ry0(-pi/4),H0,CX10,H0,Ry0(pi/4)` (1 CX).
5. Final two reflections:
   `H2,CX12,H2,Ry2(3*pi/8),CX32,Ry2(-3*pi/8)` (2 CX).

T0 is exact phase(pi/4), represented by U3(0,0,pi/4). The N3 target Rz
conjugation retains controlled(-X) exactly; its active minus sign matters.
The entire list has30 gates and11 CNOTs. Its CNOT subsequence is

`02,13,01,23,02,32,02,32,10,12,32`.

This exactly matches run409 deletion ordinal9, source position29, removing
the first CX32 of the previous final RCCX selector. The endpoint's numeric
angles were not rounded or used to derive these rational-pi angles.

## Predicted diagonal

The u0v0 branch now retains CZ instead of the old final-Y-permuted diagonal;
the remaining branch roots are unchanged:

- uv=00: `(-1)^(x*y)`.
- uv=10: `(-i)^y (-1)^x`.
- uv=01: `i^x (-1)^y`.
- uv=11: `f_i(y) (-1)^x`.

In q0q1q2q3 order, the predicted roots are

`[1,1,1,-1,1,i,-i,1,1,i,-1,-i,-1,-i,i,-1]`.

Multiplicities are6,4,3,3 for1,-1,i,-i. All angles are rational pi, with
pi/8 introducing conductor32 for a direct cyclotomic checker.

## Review and scope

Independent verifier and coordinator were sent the reflection product,
conditional-basis identity, exact topology, and predicted diagonal before
any candidate matrix evaluation. Adversarial review checked L's Bloch signs,
the N4prime half-angle, B's CX conjugation, chronological B Z versus Z B,
N3's active sign, all four inactive/active branches and honest11-CX cost.
No blocking issue found by hand. Full gate-list and independent exact
arithmetic checks remain necessary; near-hit409 alone is not certification.
The independent verifier subsequently confirmed N4prime, both final
reflections, all four complete sector operators including the d11 scalar,
the predicted diagonal, and conductor32 by hand. The coordinator's source
candidate pack is to bind these identities before any reserved computation.
