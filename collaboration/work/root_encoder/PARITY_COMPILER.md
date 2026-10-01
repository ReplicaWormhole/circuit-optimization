# Shared-parity compiler for a conditioned phase

Board49, constructive derivation; no optimizer or candidate certification.
For m chosen bits and binary word b, Pi_b=product_j(I+(-1)^b_j Z_j)/2.
Thus exp(i theta Pi_b) equals, up to the scalar exp(i theta/2^m), the product
over nonempty subsets S of exp(i theta (-1)^sum_S b_j Z_S /2^m).
Implement each factor by Rz(phi_S) on a target storing parity S, with
phi_S=-theta (-1)^sum_S b_j /2^(m-1). All factors commute.

Order selected wires, and group each subset by its highest-index wire q.
For this group take wire q as target, enumerate lower-wire subsets in cyclic
binary-reflected Gray order, and apply the corresponding Rz on its accumulated
parity. Each Gray transition flips one lower control into target via CX. Close
the cycle to restore target at the group's end; no work bit or ancilla is used.
For q>=1 this uses2^q CX. The q=0 group is one local Rz and no CX. Therefore
total cost sum_{q=1}^{m-1}2^q=2^m-2. Every nonempty subset appears once and
each parity ladder closes, proving the phase polynomial above on every basis
state. Negative controls change signs of angles only.

For theta=pi and all controls1 plus target1, conjugating this phase by H on
the target gives controlledX, up to an irrelevant input-independent scalar.
Thus CCX uses6 CX, CCCX uses14 CX, arbitrary positive/negative controls
included. This is an explicit compilation upper bound, not a minimality
claim. A controlledH follows by target-local conjugation: H=Ry(pi/4) Z
Ry(-pi/4), so k-controlledH has the same2^(k+1)-2 CX upper bound.

For the orbit construction, conditionalF4 under two sector bits requires two
two-controlledH gates(6 each) and a four-bit projector phase(theta=pi/2,
14CX), total26. Pair mixing uses three-controlledH14. Three separate F4
sectors cost78 plus14=92, improving only the previous loose268 block bound.
Alternatively apply global bit-reversedF4(2CX), correct it conditionally in
the00sector with its inverse(26), then apply the pairH(14), total42. Output
bit reversal is allowed within each sector because output eigenvalue order is
free; the inverse correction must use the SAME bit-reversedF4 as applied
globally. Generic permutation bound becomes15*7*14=1470, hence the complete
abstract family has conservative1512CX compiled bound with this correction.
Still far above13; no smaller circuit or new incumbent is claimed.

Independent reviewer must check phase signs, closed Gray cycle, global scalar,
controlledH basis and the mixed-sector composition before these compiler costs
are used in search ranking. Native multi-control gates are not free CX gates.

Independent reviewer checked these points and implemented the compiler in its
owned directory. Focused exact Fraction phase-polynomial tests cover emitted
projector identities. The native67CX full circuit (P25 plus mixing42) passes
separately reserved independent cyclotomic matrix certificate382, conductor64,
and numerical checking. This validates the construction, not a competitive
incumbent or a minimum-cost compiler. The old381 catalog remains unchanged.
