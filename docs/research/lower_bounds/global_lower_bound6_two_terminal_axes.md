# Two terminal opposite axes in a five-CNOT schedule

**Historical intermediate screen.** `global_lower_bound6_complete_filter.md`
now excludes every five-CNOT schedule and establishes the global interval
\(6\leq C_{\min}\leq14\). The counts below describe this earlier stage.

The global bound for an ancilla-free exact diagonalizer of the four-qubit
shift remains (5\leq C_{\min}\leq14). This note strengthens the necessary
schedule filter of `global_lower_bound6_terminal_filter.md`. It makes no
assumption about one-qubit gates, CNOT directions, connectivity, output
eigenvalue order, or choices inside degenerate eigenspaces.

Write (U^\dagger V_4U=D), with (D) diagonal. A wire (j) is *selected*
when it has degree two and its second incident CNOT is on an opposite
physical pair (e=\{j,k\}), with no later CNOT touching either endpoint
of (e). The first incident CNOT has a commuting input Pauli axis (p_j).
That axis passes through the first CNOT and all intervening CNOTs on other
wires, so (A_j=Up_jU^\dagger) is a one-CNOT image supported on (e).
The exact opposite-pair spectral Gram classification in
`global_lower_bound5_spectral_gram.md` says that every nonzero such image
satisfying the spectral-path identity is proportional to
\(\sigma(\mathbf n)_j\sigma(\mathbf n)_k\). Consequently
\([A_j,V_4^2]=0\), since the half-shift exchanges the two opposite wires.

The selected (p_j) has a nonzero off-diagonal entry. If it were diagonal
on input wire (j), it would equal (\pm Z_j), and (A_j) would commute
with (V_4). A non-scalar operator commuting with the four-cycle cannot
have proper physical support (Lemma 10 of the supplied paper); (A_j) is
supported on only two wires. Therefore
\([p_j,D^2]=U^\dagger[A_j,V_4^2]U=0\) and, because (D^2) is diagonal,
its eigenvalue is constant under flipping input bit (j).

If **two distinct degree-two wires** (j) and (\ell) are selected, the
eigenvalue of (D^2) is constant on every four-element square obtained by
flipping bits (j) and (\ell). Each multiplicity of (D^2) must be a
multiple of four. But (V_4^2) swaps (0\leftrightarrow2) and
(1\leftrightarrow3), so (\operatorname{tr}V_4^2=4) and its (+1,-1)
multiplicities are (10,6). This contradiction excludes the schedule.
The two selected wires may belong to **different** terminal opposite edges;
that is the new case beyond run 190.

`global_lower_bound6_two_terminal_axes.py` checks this combinatorial
criterion on all 7,776 ordered unoriented pair schedules, after the exact
necessary conditions from run 190. Among its 344 survivors, the new lemma
excludes 88. The remaining 256 are:

| Degree sequence and edge multiplicities | Remaining schedules |
|---|---:|
| (3322), five simple edges | 128 |
| (3322), one doubled edge and three simple edges | 88 |
| (4222), one doubled edge and three simple edges | 40 |

The (3322) graph with two doubled edges and one bridge is fully excluded.
These counts concern chronological unoriented pair schedules, not locally
inequivalent circuits. The result is a rigorous subclass exclusion, not a
six-CNOT lower bound. A useful next invariant is the spectral-path condition
for a selected degree-two axis whose second CNOT is followed by another
CNOT on its partner: its image can then have three-wire support, beyond the
opposite-pair Gram classification.
