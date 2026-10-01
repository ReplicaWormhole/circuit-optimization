# No Clifford-only cycle diagonalizer

Board43, structural derivation only. Qubit0 MSB, V right four-cycle;
if U V U†=D diagonal, U cannot be Clifford. This restricts gate content,
not CNOT count. Arbitrary one-qubit rotations remain allowed in the search.

Use the binary symplectic action of a Clifford on Pauli labels modulo phase,
in coordinates(x,z) in F2^4 direct-sum F2^4. V permutes wires by a four-cycle
P, so M_V=diag(P,P). In characteristic2,
(P-I)^2=P^2+I is nonzero: P² exchanges opposite wires, and P² e0+e0
equals e2+e0, not zero. Thus (M_V-I)^2 is nonzero.

Any diagonal Clifford D commutes with each Z_j, so its symplectic action
fixes every(0,e_j) individually. Writing its blocks as[A B; C E] forces
B=0,E=I. Symplectic preservation of the X/Z pairing then forces A=I;
C is symmetric. Hence M_D=[I 0; C I], giving (M_D-I)^2=0.
This includes arbitrary global phase and diagonal Pauli signs, which have
no effect on the binary action.

If U were Clifford, conjugation would give M_D=M_U M_V M_U^-1.
The square of M_D-I would be M_U(M_V-I)^2 M_U^-1 and could not vanish.
Contradiction. In a circuit built from CNOTs and local Clifford gates,
every gate and the entire U are Clifford. Thus at least one genuinely
non-Clifford local gate is required in any such gate-list diagonalizer.
This does not forbid Clifford entanglers with arbitrary local gates, does
not count non-Clifford gates optimally and supplies no new CX lower bound.

For an orbit-address construction, this explains why an affine Clifford
encoder and Clifford-only orbit mixing cannot complete the task. Non-Clifford
Fourier phases or nonlinear reversible encoding must enter somewhere. A
conditional Fourier implementation must still be costed in the requested
CX model; this obstruction alone constructs no smaller circuit.

Independent verifier reviewed the fixed-Z block constraints, symplectic
pairing, similarity invariant and SU2/global-phase scope. No defect found.
The statement is a gate-content restriction, not a lower-bound improvement.
