# Exact 17-CNOT boundary fusion

The chronological gate list is `algebraic_boundary_candidate.json`. Qubit 0 is
the most significant bit. It diagonalizes the four-qubit right cyclic shift
defined in `README.md`; it is not claimed to diagonalize total spin.

In the notation of *Exact ancilla-free diagonalization of qubit cycles*,
Appendix A, let \(P=Z_0Z_1Z_2Z_3\) and \(\alpha_j=\pi j/4\). Equation (A.2) gives
the diagonal parity correction, up to global phase, as

\[
B_4 = \exp(i3\pi P/8)
  \prod_{j=1}^3 R_{Z_j}(\alpha_j/2)
                    R_{PZ_j}(\alpha_j/2).
\]

Here \(R_M(\theta)=\exp(-i\theta M/2)\). Delete the full-parity factor to
obtain \(B'_4=B_4\exp(-i3\pi P/8)\). Since \(P\) is invariant under the cyclic
permutation, \([P,V_4]=0\), so \(F_4B'_4\) diagonalizes \(V_4\) exactly whenever
the paper's \(F_4B_4\) does.

The first Fourier gate in `baseline_18.json` is \(CZ_{1,2}\), followed by the
remaining nine-CNOT Fourier circuit \(F_{\rm tail}\). The exact identity

\[
CZ_{1,2}=e^{i\pi/4}R_{Z_1}(\pi/2)R_{Z_2}(\pi/2)
                         R_{Z_1Z_2}(-\pi/2)
\]

shows that \(CZ_{1,2}B'_4\), up to global phase, needs local \(Z\) rotations
of angles \(5\pi/8,3\pi/4,3\pi/8\) on qubits 1, 2, 3 and four nonlocal parity
rotations:

| parity mask, ordered 0-3 | angle |
| --- | ---: |
| 1011 = \(Z_0Z_2Z_3\) | \(\pi/8\) |
| 1101 = \(Z_0Z_1Z_3\) | \(\pi/4\) |
| 1110 = \(Z_0Z_1Z_2\) | \(3\pi/8\) |
| 0110 = \(Z_1Z_2\) | \(-\pi/2\) |

The following closed eight-CNOT parity network exposes every required mask.
Each row is the target-wire mask immediately after the CNOT; apply the listed
rotation there. The wire masks start at `(1000,0100,0010,0001)` and return to
that tuple at the end.

| step | CNOT | target mask | rotation |
| ---: | --- | --- | --- |
| 1 | 0→1 | 1100 | — |
| 2 | 1→3 | 1101 | \(R_z^{(3)}(\pi/4)\) |
| 3 | 2→1 | 1110 | \(R_z^{(1)}(3\pi/8)\) |
| 4 | 0→1 | 0110 | \(R_z^{(1)}(-\pi/2)\) |
| 5 | 1→3 | 1011 | \(R_z^{(3)}(\pi/8)\) |
| 6 | 0→3 | 0011 | — |
| 7 | 2→1 | 0100 | — |
| 8 | 2→3 | 0001 | — |

Conjugation by the preceding CNOTs turns each listed physical \(R_z\) into
the corresponding parity rotation. As the network restores the wire masks,
its only net action is the product of those four commuting rotations. Thus
the candidate gate list implements \(F_{\rm tail}CZ_{1,2}B'_4\) up to global
phase and has \(8+9=17\) CNOTs. This is an exact circuit identity, with no
optimized numerical angles.

`algebraic_parity_bfs.py` generated the network by breadth-first search over
closed four-wire CNOT parity bases and visited-mask flags. This search proves
only that eight CNOTs is shortest *within that parity-network formulation*;
the 17-CNOT result does not depend on its optimality. The local checker reports
maximum off-diagonal error \(4.44\times10^{-16}\). An independent Qiskit
simulation reports \(3.93\times10^{-16}\); its transpilation retains 17 CNOTs.
The file `algebraic_qiskit_result.json` has the complete numerical output.
