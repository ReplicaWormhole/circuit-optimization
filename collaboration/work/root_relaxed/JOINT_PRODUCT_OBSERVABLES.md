# Joint product-observable constraint (structural claim for review)

Let A_j=U† Z_j U, j=0,1,2,3, for a four-qubit diagonalizer UV=DU.
The A_j are commuting, traceless Hermitian involutions, commute with V,
and all nonempty products have trace zero. Thus their four generators are
independent (no nonempty product equals either signed identity).

Suppose each A_j has operator Schmidt rank one across every single-wire
cut. Successive factorization then gives a product of four nonzero local
operators. Hermiticity permits Hermitian factors after scalar rescaling;
A_j²=I normalizes each factor to a Hermitian involution. Uniqueness of a
nonzero simple tensor and cyclic invariance make the four local factors
proportional. Consequently A_j=±Q_j⊗Q_j⊗Q_j⊗Q_j. Tracelessness excludes
Q_j=±I, so Q_j=n_j·sigma for a real unit vector n_j.

For two such axes separated by angle theta in[0,pi], QR has eigenvalues
exp(±i theta). The global products commute only if (QR)^tensor4 equals
(RQ)^tensor4. A tensor eigenvector with three plus signs and one minus
sign has phase exp(2i theta); equality with its inverse requires
exp(4i theta)=1. Thus axes are parallel/antiparallel or orthogonal.
These conditions are also sufficient. At most three distinct orthogonal
axis directions occur; their global four-fold Pauli strings obey
Q1^tensor4 Q2^tensor4=Q3^tensor4 (up to any chosen generator signs),
because i^4=1. They generate at most two independent involutions.

Therefore four independent A_j cannot all have rank one across every
single-wire cut. At least two of the four output observables have rank>1 on some
single-wire cut, because any subset of these four generators is independent.
This is a joint necessity, whereas a single commuting balanced observable
can be the rank-one globalZ. It is NOT a six-CX exclusion: actual gate
reachability would have to force the conflicting rank pattern. Independent
algebraic and verifier reviews concur on the argument and its restricted
scope; direct exchanges and board handoffs record that review.
