# Necessary-condition screen for six-CNOT shift diagonalizers

Let \(V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle\), and let an
ancilla-free circuit \(U\) with arbitrary local one-qubit unitaries and
CNOTs on arbitrary pairs obey \(U^\dagger V_4U=D\), where \(D\) is
diagonal. This note applies the proven obstructions of
`global_lower_bound6_complete_filter.md` to **all** six-CNOT pair
schedules. It is a necessary-condition screen, **not a seven-CNOT lower
bound**. It makes no restriction on CNOT directions, output eigenvalue
order, or bases within degenerate eigenspaces.

## Tests and why they are necessary

Enumerate the \(6^6=46,656\) chronological sequences of unordered
physical CNOT pairs. Every input wire needs at least two CNOT incidences
by `global_lower_bound4_spectral_path.md`. Also the forward light cone
of each input wire must reach all four wires: \(UZ_jU^\dagger\) commutes
with \(V_4\), and every non-scalar operator commuting with the cyclic
shift has full physical support.

For a degree-two input wire \(j\), choose a Pauli input axis that
commutes through its first incident CNOT. Immediately before the
second incident CNOT, its image is one-body on \(j\). Let the second
CNOT act on \(j,k\). Starting from potential support \(\{j,k\}\),
scan later CNOTs in order. A CNOT is *active* if it intersects the
current potential support; add both of its endpoints to the support.
A disjoint CNOT acts trivially on the selected operator. Local gates
do not enlarge support.

- With zero active later CNOTs, the selected final image is a
  one-CNOT image on \(j,k\). The adjacent-pair spectral Gram theorem
  excludes adjacent \(j,k\); the half-shift trace identity excludes
  opposite \(j,k\). Both are proved in the six-CNOT lower-bound note.
- With exactly one active later CNOT, it must act on \(k,\ell\): wire
  \(j\) has degree two. Disjoint CNOTs, wherever interleaved, act
  trivially on this selected image. Its final image is therefore in
  the same 21-dimensional two-CNOT-chain overapproximation used by
  the three certificates referenced in
  `global_lower_bound6_complete_filter.md`. The two adjacent-first
  physical geometries are excluded by their exact rank-one spectral
  certificates; opposite-first is excluded by the exact rank-one
  classification, involution identity, and half-shift trace identity.

The screen applies a test only when its antecedent holds. For each
schedule excluded by a tail test, one degree-two selected axis yields
a contradiction for **every** choice of local gates and CNOT directions.

## Exact enumeration

`global_lower_bound7_schedule_filter.py` performs the complete finite
enumeration and writes
`global_lower_bound7_schedule_filter_result.json` (ledger run 224).
The first applicable test sets the exclusive count:

| First applicable test | Six-CNOT schedules |
|---|---:|
| Some wire has degree below two | 19,506 |
| Some forward light cone is proper | 12,198 |
| Degree-two selected axis has terminal adjacent pair | 5,064 |
| Degree-two selected axis has terminal opposite pair | 2,532 |
| One active later CNOT, adjacent chain with opposite endpoints | 1,260 |
| One active later CNOT, adjacent V geometry | 1,260 |
| One active later CNOT, opposite-first geometry | 1,260 |
| **Survive these tests** | **3,576** |
| **Total** | **46,656** |

The survivor degree sequences, written in decreasing order, are:

| Degree sequence | Schedules |
|---|---:|
| 3333 | 1,416 |
| 4332 | 1,896 |
| 4422 | 144 |
| 5322 | 120 |

The 1,416 degree-3333 schedules have no degree-two wire, so this
selected-axis tail screen cannot address them. For the remaining
survivors, the counts of active later CNOTs on their degree-two
wires are recorded in the result JSON. The script also reruns its
screen on all 7,776 five-CNOT schedules and asserts that none survive.
This regression checks the extension against the established global
six-CNOT lower bound; its exclusive category counts differ because the
one-active test can now apply when an irrelevant CNOT is interleaved.

Survival does not imply that any six-CNOT circuit diagonalizes \(V_4\).
The established global interval remains \(6\le C_{\min}\le14\). To
prove \(C_{\min}\ge7\), all 3,576 surviving pair schedules must be
excluded by tests valid for arbitrary local gates and both CNOT
directions, or by another argument covering them all.
