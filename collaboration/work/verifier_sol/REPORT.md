# Independent static selector audit

Board121; verifier_sol. Artifact types: math, computation provenance and workflow.
All work in this report is read-only source/log inspection and by-hand algebra.
No circuit matrix, objective, search, optimizer or exact-check invocation was
performed by this verifier. Root owns reservations and shared commits.

## Saved run405 comparator

The saved complete candidate reconstructs exactly from the frozen source's
literal gate builders, with50 chronological gates:19CX,17U3,8H,4Ry,2X.
The decomposition is Bellprefix2 + router2 + CSdag2 + CCX6 + CH1 + CCH6.
The saved file, config and script hashes match result.json:

- script: a6ae1b4cbd754e6f35f308ed667034011e9fb5eefc9ab89999d5edc433247bb9
- config: f42ac287926163c2f61058f9aa46af6ea3bccdf06eace6a812126d8fc29c0d93
- saved candidate: 9b299133bd7979036beaf696950f8a1b1e5959fb4d3b8507a30af9f81a84c0c5
- ledger canonical candidate: 7dfbeabdee429c42a1bc4130ceee619402ec0c7f080abc68bd639c5c053fd088

The last two byte hashes differ only by JSON formatting: parsed objects,
including metadata and the entire chronological gate word, are equal.
Current native_gate_check.py/check_circuit.py/exact_check.py hashes equal the
frozen config. Ledger405 is complete and its code/config/candidate binding
agrees. The saved wrapper/checker stderr files are empty; both returncodes
are0 and neither checker timed out.

The numerical log reports maxoff3.1060360143791505e-16 and unitarity error
1.7763568394002505e-15. The exact log reports conductor16 and16 row labels,
with counts(6,3,4,3) in the exact_check.py root order(1,i,-1,-i). Inspection of
exact_check.py confirms that its certificate checks each row's full16-entry
right-shift eigenvector relation in chronological/MSB-zero conventions.
This independently audits the saved exact certificate's input/provenance;
it does not rerun or replace that certificate. The comparator exceeds12CX
and does not improve the accepted13CX incumbent.

The frozen CH source comment reverses the rotation signs, as already noted
in run405/PLAN.md. Its actual chronological Ry(-pi/4),H,CX,H,Ry(+pi/4)
word correctly gives active Ry(+pi/4) Z Ry(-pi/4)=H. The saved gate word,
not that comment, is authoritative.

## Independent by-hand counterexample for the literal12CX variants

This confirms the sector argument supplied directly by algebraic_sol and
applies only to the literal RCCX/CZ substitution choices described below.
Use A=(q0,q2)=(a0,a1), d=(q1,q3)=(d0,d1), MSB qubit0 and chronological
gates. The Bell decoder and router give T_d=CZ_A X_A^d. The proposed word
uses Btilde=CZ_A, P=RCCX(q3,q0->q2) with optional control swap/adjoint,
CH(q1->q0), then a negative-q1-control CCH made by target conjugating an
RCCX(q1,q3->q2), again with optional control swap/adjoint.

The literal chronological H,T,CX_b,Tdag,CX_a,T,CX_b,Tdag,H has target
blocks I,I,Z,Y for controls(a,b)=00,01,10,11, respectively. These blocks
follow directly by tracking a computational target bit t between the two H
gates: the net flip is a and the T exponent is
t-(t xor b)+(t xor a xor b)-(t xor a). This gives inner I,I,X,-Y
and hence outer I,I,Z,Y. These blocks are Hermitian, so the full RCCX is
Hermitian/involutory and its adjoint
choice changes no unitary. Four control-order pairs exhaust the distinct
unitaries among the16 declared choices.

Fix d=(1,0). In the final CCH, the surrounding X(q1) changes q1 to0;
therefore both controls are0 for either control order. Its RCCX is I and
its target conjugations cancel exactly. In P, order(a,b)=(q3,q0) gives I;
order(a,b)=(q0,q3) gives Z_q2 when a0=1, namely CZ_A. Consequently the
remaining branch basis change is C_A=H_q0 CZ_A or H_q0.

Restrict further to a1=1. The target is T10=Z_q0 X_q0=iY_q0. Conjugating
by Z or H only flips the sign of Y, since ZYZ=-Y and HYH=-Y. Therefore
C_A T10 C_A^dagger is plus or minus iY_q0 in either case, with
computational off-diagonal entries of magnitude1. It is not diagonal.

This exact counterexample invalidates every literal variant irrespective
of forward/adjoint bits. It does not exclude other phase corrections,
other relative-phase selector decompositions, other12CX circuits or the
full98 family. No numerical screen is necessary to establish this narrow
failure, although a separately reserved finite screen could corroborate it.

## Review and cleanup

Adversarial review used ~/.codex/adversarial_review.md. The review checked
conjugation order, negative control behavior, both P control orders, the
exact root-order convention, hashes versus semantic JSON identity, and
claim strength. No unresolved defect was found in this report. A fresh
numerical_sol script/config is still pending static readiness review;
no prior unfinished script is approved by this report. This directory is
intentional research evidence; other dirty/untracked work was preserved.
