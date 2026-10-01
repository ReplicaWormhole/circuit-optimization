# An exact selected-axis witness for the adjacent-centered tail

This result concerns **one selected input Pauli axis** in the
adjacent-then-adjacent degree-three two-active-CNOT tail. It does not
construct a diagonalizer of the four-qubit shift and does not exclude
any six-CNOT schedule.

Let \(j=0,k=1,\ell=3\), and choose

\[
p=(X_0-Y_0+Z_0)/\sqrt3,\qquad
W=\mathrm{CZ}_{03}\,H_0\,\mathrm{CNOT}_{0\to1}.
\]

The controlled-Z uses one CNOT and local Hadamards, so \(W\) has two
entangling gates on the ordered pairs \((0,1),(0,3)\). Direct Pauli
conjugation gives

\[
A=WpW^\dagger
 = (Z_0X_1+X_0Z_3+Y_0X_1Z_3)/\sqrt3.
\]

The three displayed Pauli strings anticommute pairwise, proving
\(A^2=I\). Exact 16-by-16 SymPy multiplication in
`global_lower_bound7_degree3_centered_aa_witness.py` also verifies
that \(P_a A P_b A P_c=0\) for every triple of distinct shift
eigenvalue labels. Thus \(A\) satisfies the single-axis spectral-path
identity and is an actual two-CNOT image, not merely a vector in
the 24-dimensional overapproximation.

The 14-real-variable common-first-partner-axis system in
`global_lower_bound7_degree3_centered_aa_common_axis.py` includes
the spectral-path Gram, \(A^2=I\), and unit normalization of the
first partner's Bloch axis. Its exact rational Gröbner basis has
real solutions, including this witness (ledger run 248). The exact
circuit and spectral check is ledger run 253; failed run 252 used
unsimplified SymPy zero comparison and was superseded.

The exact shift trace moments are

\[
\bigl(\operatorname{tr}(V_4A),\operatorname{tr}(V_4^2A),
\operatorname{tr}(V_4^3A)\bigr)
 = (2i/\sqrt3,0,-2i/\sqrt3).
\]

For this input axis, \(p_z=1/\sqrt3\). The finite eigenvalue-count
test in `global_lower_bound7_degree3_centered_aa_trace_assignment.py`
finds a bit-\(0\) assignment with counts \((3,2,2,1)\) for
eigenvalues \((1,i,-1,-i)\), and bit-\(1\) counts \((3,1,2,2)\).
These use the correct total multiplicities \((6,3,4,3)\) and match
all three trace moments exactly (ledger run 256). This is only a
necessary one-bit consistency check, not a full diagonal eigenvalue
assignment or a circuit construction.

Consequently the combined **single-axis** tests used here—spectral
path, Hermitian involution, common first-partner axis, and these
trace moments—cannot exclude the entire adjacent-centered physical
shape. Constraints involving several input axes or all diagonal
Pauli masks remain possible. The sound six-CNOT schedule screen
still leaves **3,208** ordered unoriented pair schedules.
