# Independent four-CX merge check

Read-only by-hand check of algebraic_sol's proposed merge; no matrices,
search or exact-check invocation. The repaired first P must use controls
(x=q0,v=q3), target y=q2, with chronological S_y,RCCX,Sdag_y.
Let B=CSdag on (x,y). PB has target blocks I,I,S,XSdag in
(x,v)=00,01,10,11. P has I for x=0, Z for x1v0 and X for x1v1;
B has I for x=0 and Sdag for x=1. Thus Z Sdag=S proves the block table.

Strip the output T_x=diag(1,exp(i*pi/4)). The remaining x1v0 block is
A=exp(-i*pi/4)S=Rz(pi/2), and x1v1 block is
C=exp(-i*pi/4)X Sdag=-i(X+Y)/sqrt2.

Set N1=(-X+Y)/sqrt2, N2=N4=(-X+Z)/sqrt2, N3=-X. The chronological
controlled-N sequence N1(x),N2(v),N3(x),N4(v) gives blocks
I,N4N2,N3N1,N4N3N2N1. Direct Pauli products give:

- N4N2=I;
- N3N1=(I-iZ)/sqrt2=Rz(pi/2)=A;
- N2(-X)N2=Z, so N4N3N2N1=ZN1=-i(X+Y)/sqrt2=C.

Appending T_x therefore reproduces PB exactly, including the relative
phases between control sectors. Each controlled N is one CX with target
conjugations. One complete chronological block is:

1. Rz_y(-3pi/4), CX(x->y), Rz_y(3pi/4).
2. Ry_y(3pi/4), CX(v->y), Ry_y(-3pi/4).
3. Rz_y(-pi), CX(x->y), Rz_y(pi).
4. Ry_y(3pi/4), CX(v->y), Ry_y(-3pi/4).
5. T_x as exact U3(0,0,pi/4).

N1=Rz(3pi/4) X Rz(-3pi/4), N2=Ry(-3pi/4) X Ry(3pi/4),
and N3=Rz(pi) X Rz(-pi)=-X, confirming the chronological signs.
The third block can equivalently be Z_x CX(x->y), but the target
conjugation form matches the uniform elementary construction.

Together with Bellprefix2CX, router2CX, CH1CX and repaired final CCH3CX,
the cost is2+2+4+1+3=12CX. The repaired sector derivation applies because
the four-CX block is exactly PB for the required swapped firstP order.
This is a structural12CX proposal until a complete saved gate list has
independent validation and exact certification. No incumbent has changed.

Adversarial review checked all four block products, T's relative phase,
chronological signs and the firstP control-order requirement. No unresolved
issue found in this by-hand derivation; full saved-list certification pending.
