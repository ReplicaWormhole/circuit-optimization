# Exact absorption of control-wire interior local gates

Prefix chronology L0 CX01 A CX21 B CX31 L1, with four-wire arbitrary
product-local L0,A,B,L1. Q0 is most significant, no ancillas; directed CX.
The three CX controls are distinct and have common targetq1.

A gate on any controlwire can be moved to a boundary without crossing its
own controlled gate. In detail:
- A(q0) crosses CX21,CX31 forward and absorbs into L1(q0); combine B(q0).
- A(q2) crosses CX01 backward and absorbs into L0(q2).
- A(q3) crosses CX01 backward; B(q3) crosses CX21,CX01 backward, combining
  with A(q3) on that wire before absorption into L0(q3).
- B(q2) crosses CX31 forward and absorbs into L1(q2).
All other local factors are on disjoint wires. These moves preserve exact
unitaries, including local phases, by chronological multiplication.
Thus the family can be represented with only A(q1),B(q1) in the two interior
layers, six rather than24 additional coordinates. This is an equality of
realizable families with free boundary layers, not permission to freeze
boundary angles while deleting optimizer coordinates.

Further Rx rotations on targetq1 commute through each CX. Their Euler/gauge
freedom creates redundancies; no dimension or Hessian nullity is asserted
without tracking the corresponding changes in boundary coordinates. A
projected24-coordinate Hessian with boundaries fixed can have nonzero modes
along directions that are redundant in the full parameterization. Negative
curvature still yields a valid finite-dimensional descent direction; a
nonnegative projected Hessian certifies neither full116 PSD nor a minimum.

Algebraic and verifier independent reviews confirm the commutation proof.
An equivalent future family has98 coordinates:84 boundary/exterior local
coordinates, six target-interior coordinates and eight F angles. This count
is a parameterization ceiling, not a dimension or optimal-gate proof. No numerical experiment or certificate
executed by this structural note.
