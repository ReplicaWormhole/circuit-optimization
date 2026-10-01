# Full support for halfshift-commuting input-axis images

Structural consequence, board32. No new exact search or schedule enumeration.
Let V be the four-wire right cycle, let U V U†=D diagonal, and let
A=U† p_j U for a unit Pauli axis p on one output wire j. This is the inverse-U
translation of the convention in the existing opposite-pair theorem, which
uses U_old† V U_old=D and A=U_old p_j U_old†. Set U_old=U†.

Define the minimal support S(A) as the union of wires occurring nontrivially
in the Pauli expansion of A. Uniqueness of the Pauli basis implies that
conjugating by a wire permutation sends S(A) to its permuted set. If
[A,V²]=0, S(A) is invariant under the permutation (0 2)(1 3). Its invariant
subsets are empty, {0,2}, {1,3}, and all four wires. Empty support means A is
scalar; this is impossible for a traceless Hermitian involution obtained by
conjugating a unit Pauli axis. The existing exact opposite-pair input-axis
obstruction forbids both two-wire alternatives. Therefore S(A) is all four.

Consequently no proper-support input-axis image can commute with the halfshift.
This is stronger in support scope than a trace-zero statement, but is already
a consequence of the proved opposite-pair obstruction, not a new independent
lower-bound theorem. If a schedule-specific argument forces a selected axis
to commute with V² while its support cone omits a wire, the schedule is
excluded. No such additional forcing argument is established here and the
3024 six-CX necessary-screen survivors remain unchanged.

Assumptions are essential: A must be a one-wire Pauli-axis image of a genuine
diagonalizer. Arbitrary operators commuting with V² can have opposite-pair
support (for example SWAP02), so the conclusion is not a commutant classification.
The theorem imported is `global_lower_bound7_opposite_pair_support_certificate.md`,
with exact evidence in completed runs356–357. Its historical upper bound14 is
superseded by the accepted exact13 circuit; its opposite-pair proof is unchanged.

Independent verifier reviewed the convention translation, support alternatives
and applicability of the existing opposite-pair theorem directly. No defect
found. Neither the old certificate nor a new symbolic search was executed.
