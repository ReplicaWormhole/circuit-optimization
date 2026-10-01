# Independent exact local-rewrite review

Hypothesis52; no optimizer or parent382 recertification. Exact67sourcehash
d595df83cc218ab4a532cdc2c51ab90a19dd68503fccc8d10311a23f0608ed6d matches
completed382 config. Actual frozen rewrite rules independently reviewed
before383reservation/execution: disjoint gates, sameaxis rotations, controlRz,
targetRx and CX-CX c1!=t2,c2!=t1 commute; identicalH/CX cancel; rational-pi
sameaxis angles add exactly without modulo2pi global-phase dropping.

Independent replay of all5 actual trace rules verifies each crossed gate,
source operands/replacement and final native list exactly. CandidateSHA
99742270d05745155e74379602293bfb7bf2c19782d57def349521a570a000f5,
149→143gates,67CX unchanged. Independently assembled NumPy matrices differ
from source at most8.40036e-16 without phase alignment; candidate maxoff
7.62780e-16, unitarity2.44e-15. Qubit0 MSB chronological UV=DU. Passing
numeric evidence and exact rewrite identities do not constitute a fresh
full UV certificate; root notified to reserve one new exact check separately.

Algebraic orbit-phase criterion independently reviewed: for fixedorbit basis,
weighted successor coefficients r_(j+1)conj(r_j) are constantρ iff r_j=r0ρ^j,
ρ^m=1. PositiveDFT sendsρS toρdiagω^k. Actual bit-reversedquartet root order
[0,2,1,3] shifts to[s,s+2,s+1,s+3]mod4; pair shifts[2s,2+2s]mod4;
singleton globalphase leaves label0.128discrete sectorchoices and six orbit
phaseconstants are within fixedbasis scope, not arbitrary degeneracy gauges.
A proposed12CX family count3+4*2+1 is arithmetic only until actual gates pass.

Command OPENBLAS_NUM_THREADS=1 python3
collaboration/work/verifier_label12/joint_simplify/check.py passes. Reviewed
bounded383 one120second/100sweep/50000rulecheck plan; actualfixpoint1225checks
under0.001second. Adversarial review found no confirmed defect. No CXreduction,
incumbent promotion or lowerbound claim; hybridnativeF counts separate.
All artifacts retained intentionally; coordinator owns commit, no push.
