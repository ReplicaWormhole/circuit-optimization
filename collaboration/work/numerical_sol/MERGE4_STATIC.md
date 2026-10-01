# Independent by-hand audit of the four-CX B/P merge

Board126. This is static algebra: no matrices were evaluated by a program,
no optimization was run, and no scientific ledger attempt was needed.
Algebraic_sol proposed the identity; numerical_sol independently checked its
four conditional target blocks and primitive gate chronology.

Let x=q0, v=q3 and y=q2. The corrected first P is chronological
S_y,RCCX(x,v;y),Sdg_y, so its target blocks in xv order00,01,10,11 are
I,I,Z,X. For B=CSdg(x,y), the chronological merged block P B has blocks
I,I,Z Sdg,X Sdg = I,I,S,X Sdg. Controls x and v are unchanged throughout.

Define unit-length Pauli reflections

- N1=(-X+Y)/sqrt2;
- N2=N4=(-X+Z)/sqrt2;
- N3=-X.

Use the four chronological controlled reflections N1 conditioned on x,
N2 conditioned on v, N3 conditioned on x, N4 conditioned on v, then T_x.
The target products before output T_x are:

| xv | Chronological activated gates | Matrix product |
| --- | --- | --- |
| 00 | none | I |
| 01 | N2,N4 | N4 N2 = I |
| 10 | N1,N3 | N3 N1 = (I-iZ)/sqrt2 = Rz(pi/2) |
| 11 | N1,N2,N3,N4 | N4 N3 N2 N1 = -i(X+Y)/sqrt2 |

For the last identity, N2(-X)N2=Z, and
Z(-X+Y)/sqrt2=(-iY-iX)/sqrt2. Multiplying the x=1 branches by the
output control phase T_x=e^(i pi/4) gives

- e^(i pi/4) Rz(pi/2) = diag(1,i)=S;
- e^(i pi/4)[-i(X+Y)/sqrt2] = X Sdg.

The x=0 branches remain I. Thus all four target blocks agree exactly with
P B, including phases, rather than merely agreeing as classical permutations.
This merged block uses four CX instead of five.

## Literal one-CX implementations

Each reflection uses a target conjugation of X, which costs one controlled-X.
For standard Rz(theta)=exp(-i theta Z/2) and Ry(theta)=exp(-i theta Y/2),
the exact chronological lists are:

1. N1 conditioned on x: Rz_y(-3pi/4),CX(x,y),Rz_y(3pi/4).
2. N2 conditioned on v: Ry_y(3pi/4),CX(v,y),Ry_y(-3pi/4).
3. N3 conditioned on x: Rz_y(-pi),CX(x,y),Rz_y(pi).
4. N4 conditioned on v: Ry_y(3pi/4),CX(v,y),Ry_y(-3pi/4).
5. Output T_x: U3_x(theta=0*pi,phi=0*pi,lam=pi/4).

Inactive-control branches are identity because each target conjugation cancels
exactly. No unrecorded local phase or uncontrolled gate is omitted.

Replacing B/P in the repaired construction gives nominal total
2(Bell)+2(router)+4(merged B/P)+1(CH)+3(final relative CCH)=12CX.
The remaining gates must use the independently reviewed repaired final wrapper
X1,Ry2(3pi/4),Rx2(pi/2),RCCX(1,3;2),Rx2(-pi/2),Ry2(-3pi/4),X1.
Qubit0 is most significant and gates are chronological. The complete12CX
gate list still requires the coordinator's reserved numerical and exact
certification. This static merged-block identity alone is not an incumbent
acceptance, a global optimum proof, or a strong Schur-transform result.

Adversarial review: all four control branches and exact scalar phases were
checked by Pauli multiplication, including inactive branches. Rotation signs
were checked against Rz X Rz†=cos(theta)X+sin(theta)Y and
Ry X Ry†=cos(theta)X-sin(theta)Z. No findings in this block derivation.
Remaining risk is transcription and complete-circuit composition, assigned to
the independent verifier. No existing files or source13 baseline were modified.
