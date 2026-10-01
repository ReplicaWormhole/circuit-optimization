# Independent by-hand review of repaired13CX selectors

Board125. This review touches math and workflow. No matrices, symbolic
search, optimizer or exact checker were run. The full saved gate list and
exact certificate remain required before treating this as a certified
construction. This is a13CX anchor, with no improvement to the accepted13CX
incumbent and no12CX result.

Use A=(x,y)=(q0,q2), d=(u,v)=(q1,q3), MSB qubit0 and chronological gates.
The Bell prefix costs2CX and router2CX. B=CSdag or CS costs2CX.
Write s=+i for CSdag or -i for CS. After B the four routed targets are
CZ for d00, s^y X_x for d10, s^x X_y for d01, and
f(x xor y) X_x X_y for d11, with f(0)=s,f(1)=1.

## First selector repair

The literal RCCX blocks are I,I,Z,Y in control order00,01,10,11.
Conjugate its target by Sdag: chronological S,RCCX,Sdag. This gives
I,I,Z,X because Sdag Y S=X and Sdag Z S=Z. It still costs3CX.
Let P be this gate with controls(v,x), or with swapped controls(x,v),
and target y.

For v=1 its target blocks indexed by x are (Z,X) or (I,X), respectively.
Both preserve x and route y to y xor x. Consequently they route
f(x xor y) to f(y). On X_x X_y their conjugations are X_x Z_y or X_x,
since the cross products K1 X K0 and K0 X K1 both equal Z or I.
For v=0 the two orders give I or CZ_A, respectively.

## Final selector repair

Keep the exact one-CX CH controlled by u on x. For the final3CX selector
use controls(c,v), where c=1-u is implemented by surrounding X(q1),
and target-conjugate positive-sign RCCX by
Vprime=Ry(-3pi/4) Rx(-pi/2).

Chronological target word is
Ry(3pi/4),Rx(pi/2),RCCX,Rx(-pi/2),Ry(-3pi/4).
The Bloch rotations show Vprime Y Vprime-dagger=H and
Vprime Z Vprime-dagger=Y: Rx(-pi/2) sends Y to -Z and Z to Y;
Ry(-3pi/4) sends -Z to (X+Z)/sqrt2 and fixes Y.
Thus the final sector actions are Y_y on d00, H_y on d01, and I on
d10,d11. The inactive error is confined to the computational CZ sector.

## Complete sector check

After P, CH, and the final selector, the targets are:

| sector | P controls(v,x) | P controls(x,v) |
|---|---|---|
| d00 | Z_x CZ_A | Z_x CZ_A |
| d10 | s^y Z_x | s^y Z_x Z_y |
| d01 | -Z_x s^x Z_y | s^x Z_y |
| d11 | f(y) Z_x Z_y | f(y) Z_x |

All entries are computational diagonal. In d00 the final Y conjugates
CZ_A to Z_x CZ_A, merely permuting its computational eigenbasis. In d01,
P controls(v,x) conjugates X_y to -Z_x X_y; the swapped order leaves X_y
unchanged. Final H_y diagonalizes either. In d10, I or CZ_A contributes
an optional Z_y while H_x diagonalizes X_x. In d11 the cross-product
calculation above and H_x give the stated diagonal forms. Both signs of B
are allowed because all s-dependent factors are already diagonal.

The complete cost is2+2+2+3+1+3=13CX. This is a by-hand diagonalization
of the stated abstract literal primitives for both B signs and both first
P control orders. The final controls must be ordered(c,v); swapping those
controls would place the inactive Y on d11 and is not covered by this table.

Adversarial review used ~/.codex/adversarial_review.md, checking chronology,
both B signs, both P orders, rotation signs, all inactive-control branches,
and the distinction between structural derivation and saved-list exact
certification. No unresolved issue found in this limited derivation.
Root must freeze/reserve any finite gate-list verification before execution.
