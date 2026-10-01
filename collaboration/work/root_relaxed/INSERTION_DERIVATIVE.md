# Local insertion derivative of the label-free loss

Assume an ancilla-free unitary U and the four-qubit right cycle V, q0 most
significant. Put B=UVU†, d_j=B_jj, and L=||offdiag(B)||F²/16.
Since V and U are unitary, L=1-sum_j |d_j|²/16.
Consider U(t)=A exp(-itP/2) C, U(0)=AC, with Hermitian single-wire Pauli P.
With K=APA†, B(t)=exp(-itK/2) B exp(itK/2), hence
B'(0)=-(i/2)[K,B], and

L'(0)=-(1/8) Re sum_j conjugate(d_j) B'(0)_jj
     =-(1/16) Im sum_j conjugate(d_j) [K,B]_jj.

This is an insertion into an actual circuit, not a free choice of arbitrary
K: A is the specified suffix. Setting two added local layers to identity
embeds an old circuit, but does not imply those24 added derivatives vanish.
Conversely, zero insertion derivatives alone neither certify a minimum nor
exclude a successful finite perturbation. A continuation that stops at an
embedded stationary point must not claim the expanded family is excluded.

The numerical researcher is asked to retain initial loss/full gradient norm
and the norm of the24 new A/B gradient components within the already reserved
fit evaluations. No additional starts or curvature experiment is authorized
by this structural note. The verifier should check the derivative sign and
scope directly; it requires no incumbent recertification.

Independent verifier review confirms the sign and actual-suffix restriction.
Run386 measured added-coordinate gradient norms6.03e-9 and2.20e-8; these
are small numerical derivatives, not exact zeros or a minimum certificate.
