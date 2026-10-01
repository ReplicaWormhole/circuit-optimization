# Independent static obstruction for the corrected13CX anchors

Board122. This is a by-hand structural review only; no matrix, symbolic
search, optimizer or exact checker was run. Numerical_sol had not frozen
config/manifest when the obstruction was found. Root owns any execution.

Scope: Bellprefix2 + R2 + B2 + P3 + CH1 + CCH3, with B=CSdag or CS;
P is the literal RCCX with either control order and either the displayed
T/Tdag signs or their complex conjugates; final CCH uses positive-sign
RCCX with either control order, negative q1 control and target conjugation
V=Ry(-pi/4)Sdag mapping Y to H. This13CX word is an anchor proposal only.

Use x=q0,y=q2,v=q3 and d=(q1,q3). Fix d=11. Before P, after either B sign,
the routed target is f(x xor y) X_x X_y where f(odd)=1 and f(even)=+i
for CSdag or -i for CS. P preserves x and maps y to y xor x. Therefore it
conjugates the parity phase into f(y).

Let epsilon=+1 for the displayed RCCX and -1 for its T-sign conjugate.
The literal RCCX target blocks are I,I,Z,epsilonY for controls00,01,10,11.
Thus for v=1 its target blocks indexed by x are:

- controls(v,x): K0=Z, K1=epsilonY;
- controls(x,v): K0=I, K1=epsilonY.

Because X_x X_y changes x, P(X_x X_y)Pdag has target operator K1 X K0
on input x=0 and K0 X K1 on input x=1. In the first control order these
are -i epsilon I and +i epsilon I; in the second they are
-i epsilon Z and +i epsilon Z. Hence the conjugated target is
-epsilon Y_x tensor f(y) in the first order and
-epsilon Y_x tensor f(y)Z_y in the second.

CH applies H_x in d=11 and HYH=-Y, leaving the same Y_x tensor factor
up to sign. The final negative-control wrapper turns q1 into0, while
q3=1. Its RCCX is I or Z_y for the two orders, so target conjugation gives
I or VZ_yVdag. This acts solely on y and cannot remove the nonzero Y_x
off-diagonal factor. The final target has form plus/minus Y_x tensor K_y
with K_y unitary, so it cannot be computational diagonal.

This independently confirms algebraic_sol's sector argument for all16
corrected13CX literal choices. There is no need for a numerical screen to
prove this restricted failure. It does not exclude different RCCXs, added
local phases, different phase-correction placements or arbitrary13/12CX
circuits. It changes no incumbent and supplies no global lower bound.

Adversarial review followed ~/.codex/adversarial_review.md: both B signs,
both P orders/T signs, target chronology, final inactive-control blocks,
and claim scope were checked. No unresolved defect found. The original
numerical_sol12CZ script remains preserved/unexecuted; no fresh13 config
or script was approved or modified by this review.
