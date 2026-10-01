# Two selected-axis test for one six-CNOT schedule

This is a bounded joint necessary-condition test for the ordered
unoriented-pair schedule

\[
(01),(23),(12),(01),(23),(03).
\]

It survives the earlier sound filters. Wire 0 has a selected axis
whose two active CNOTs are \((01),(03)\); wire 3 has one whose two
active CNOTs are \((23),(03)\). Both are adjacent-then-adjacent
centered tails. `global_lower_bound7_joint_aa_schedule_select.py`
enumerates the current 3,208 six-CNOT survivors and selects this
representative (ledger run 259).

Fix the canonical exact wire-0 image from run 253,

\[
A=(Z_0X_1+X_0Z_3+Y_0X_1Z_3)/\sqrt3.
\]

For the selected input axis on wire 3, each active CNOT preserves a
nonidentity factor on wire 3. Before using any additional collective
rotation, its final image lies in the 48-dimensional real Pauli span
with arbitrary factors on wires 0 and 2, identity on wire 1, and a
nonidentity factor on wire 3. This is an overapproximation valid for
all local gates and CNOT directions. Images of axes on different
input wires commute, so a possible second image \(B\) must obey
\([A,B]=0\).

`global_lower_bound7_joint_aa_centralizer.py` computes this linear
condition exactly in the Pauli basis (ledger run 260). The
commutator map has rank 36; its kernel inside the 48-dimensional
span has dimension **12**. A basis consists of the twelve strings
\(I_0I_1\sigma_2Z_3\) and
\(Z_0I_1\sigma_2X_3, Z_0I_1\sigma_2Y_3\), where
\(\sigma_2\in\{I,X,Y,Z\}\).

The second spectral-path and involution conditions are feasible in
this centralizer. An exact witness is

\[
B=(Z_3X_2+X_3Z_0+Y_3X_2Z_0)/\sqrt3.
\]

The three terms of each image anticommute pairwise, so
\(A^2=B^2=I\). Exact 16-by-16 SymPy arithmetic verifies
\([A,B]=0\) and all spectral-path identities
\(P_a A P_b A P_c=P_a B P_b B P_c=0\) for distinct labels.
`global_lower_bound7_joint_aa_pair_witness.py` records this check
(corrected ledger run 263). Its first version, run 261, was marked
failed because unsimplified SymPy expressions were compared to zero.
The correction simplifies the exact matrices before comparison.

This proves only that the **two selected-axis necessary conditions**
are jointly satisfiable for the fixed wire-0 witness branch. It does
not prove that one six-CNOT circuit realizes both axes, nor that the
schedule diagonalizes \(V_4\). Conversely, the tested conditions
cannot exclude this representative. Physical dihedral relabeling
transports the same feasible test to its orbit of eight schedules;
reflection conjugates \(V_4\) to \(V_4^{-1}\), which permutes the
spectral labels and preserves these identities. No statement is
made about other six-CNOT schedule classes. The global sound screen
still leaves **3,208** ordered unoriented-pair schedules.
