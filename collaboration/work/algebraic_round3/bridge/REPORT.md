# Decorated bridge derivation: hypothesis24

Structural math and workflow only: no optimizer, symbolic search or ledger run.
Let C=CNOT(0,1), E=CNOT(1,2), F=CNOT(0,2). Both interior tensor-product local layers L1,L2 must commute with C. A sufficient family is wire0 Rz, wire1 Rx, and arbitrary spectator-wire2/3 local gates. Exterior layers Lk,L3 are arbitrary. No full local-centralizer classification is claimed.

The chronological motif Lk,C,L1,E,L2,C,L3 has matrix

    L3 C L2 E L1 C Lk = L3 L2 C E C L1 Lk = L3 L2 F E L1 Lk.

Thus its chronological collapse is Lk,L1,E,F,L2,L3. In the fixed-slot template use [None,E,F], retain Lk,L1, set the middle local layer to I, and replace the next layer with L3@L2 independently per wire. Spectators need not commute with E or F; only the two C commutations are needed. SU2 products remain SU2. Euler/U3 conversion can introduce a global phase; compare complete saved matrices up to phase separately from diagonalization checks.

Canonical [01,32,12] requires two mutations to the bridge [01,12,01], followed by collapse to [12,02]. This is not an exact replacement of the incumbent. Floating commutant fits require independent numerical checks and exact angle recognition/certification before incumbent promotion.

Independent verifier confirmed both identity and fixed-slot conversion. Self adversarial review checked chronological order, spectator scope and two-mutation provenance; no confirmed issue. No smaller diagonalizer, topology exclusion or optimum claim follows. Root owns final shared review and commit.

## Follow-up: complete tensor-product local centralizer

The verifier asked whether a larger local family could preserve the bridge while retaining its CNOT cost. For exact commuting tensor-product locals the answer is no. On its two active wires,

    C = (I + Z0 + X1 - Z0 X1)/2.

Partial traces of C are I+Z on the control and I+X on the target. If (A tensor B) C (A tensor B) dagger = C, partial tracing forces A Z A dagger=Z and B X B dagger=X. Conversely those two equations imply commutation. Thus A and B are respectively Rz and Rx up to scalar phases; there are no further disconnected tensor-product components. Commutation up to a phase also gives phase1 because trace(C)=2 is nonzero. Spectator local gates remain unrestricted. Arbitrary entangling centralizers are a different family whose synthesis cost must be counted separately. This definition-level derivation was sent for independent review and does not alter the current numerical run.
