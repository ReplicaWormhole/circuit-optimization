# Canonical13 structural analysis — hypothesis17

Researcher: algebraic_round2. Coordinator persisted the returned report because
child shell execution failed. No optimizer, symbolic search or duplicate
exact13 certification was run. Findings were exchanged directly with numerical
and verifier researchers before final collection.

The thirteen directed interactions are 01,32,12,10,23,21,01,32,12,10,23,21,01.
They form two decorated six-interaction sweeps followed by01. Noncommuting local
rotations prevent applying bare-CNOT cancellation identities directly.

Write the certified circuit as U=L C01 W. Terminal locals are
L0=Ry(3pi/4), L1=Ry(3pi/4)Rz(phi), cos(phi)=-sqrt6/3, sin(phi)=sqrt3/3.
Thus A=L0 Z L0*= (X-Z)/sqrt2, B=L1 X L1*=(X+Y+Z)/sqrt3.
Deleting only terminal C01 gives U'=KU, with
K=L C01 L*= (I + A tensor I + I tensor B - A tensor B)/2.
Its column0 has K[8,0]=(1-1/sqrt3)/(2sqrt2)>0 and
K[0,0]=(1-1/sqrt2+1/sqrt3+1/sqrt6)/2>0. The recorded output labels for
rows0 and8 are -1 and+1. Since K is Hermitian, K[0,8]=K[8,0]>0.
If E=KDK* were diagonal, EK=KD would imply E[0,0]=D[0,0]=-1
from K[0,0]>0 and E[0,0]=D[8,8]=+1 from K[0,8]>0, a contradiction.
The shared row, rather than column support against the original labels, proves
the obstruction even when E may reorder labels. Coordinator corrected this
logical detail during adversarial review.
Consequently unchanged-angle terminal deletion fails, even with free output
label order. This excludes only that fixed deletion, not refitted twelve-CNOT
circuits. Coordinator independently checked the terminal gate metadata and signs.

The chronological identity C01 C12 C01 = C12 C02 follows from binary wire
updates and saves one CNOT when intervening local gates can be moved out or
arranged into compatible gadgets. Current decorated blocks do not establish
such a reduction. No smaller construction or topology exclusion was found.

Next hypothesis: jointly replace a decorated three-wire segment; freeze two
changed twelve-CNOT schedules, fit full local layers retaining eigenspace
freedom. Suggested future budget one start per schedule,300iterations and
120seconds each, one thread. This recommendation is not dispatched.
