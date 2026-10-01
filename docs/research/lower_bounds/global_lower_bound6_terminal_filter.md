# Five-CNOT topology audit after the spectral-path bound

**Historical intermediate screen.** The later trace obstruction and complete
schedule proof in `global_lower_bound6_complete_filter.md` exclude all
five-CNOT schedules. Counts below describe this earlier stage only.

This is a bounded necessary-condition audit, **not** a six-CNOT lower bound.
The global result remains \(5\leq C_{\min}\leq14\). The proof of the
five-CNOT lower bound is in global_lower_bound5_spectral_gram.md.
We allow arbitrary one-qubit gates, CNOT directions, CNOT pairs, eigenvalue
orderings, and bases in degenerate eigenspaces.

## Degree and multigraph classification

Every wire has at least two CNOT incidences by the spectral-path argument
in global_lower_bound4_spectral_path.md. Five CNOTs give ten incidences,
so the sorted degree sequence is either \((4,2,2,2)\) or \((3,3,2,2)\).
The interaction multigraph must be connected by Section 7, Lemma 10 and
the proof of Theorem 2 in the supplied *Exact ancilla-free
diagonalization of qubit cycles* paper.

For \((4,2,2,2)\), the unique connected multigraph type is a triangle
through the degree-four vertex, with a doubled edge from that vertex to
the fourth vertex. To see this, four edges touch the degree-four vertex;
the remaining edge joins two degree-two vertices. They each have one
remaining incidence to the high vertex, while the third degree-two
vertex has two parallel incidences to it.

For \((3,3,2,2)\), let \(x\) be the edge multiplicity between the two
degree-three vertices and \(y\) between the two degree-two vertices.
The cross-edge count implies \(x=y+1\). The connected possibilities are:

1. \(x=1,y=0\), with all four cross edges distinct: \(K_4\) minus the
   edge between the degree-two vertices.
2. \(x=1,y=0\), with two doubled cross edges, one from each high vertex
   to a different low vertex, plus the high-high edge.
3. \(x=2,y=1\), with a cross matching joining the two high and two low
   vertices.

The remaining \(x=3,y=2\) possibility is disconnected (triple high-high
and double low-low edges), so it cannot diagonalize \(V_4\).

## Terminal-edge obstruction

Call an edge *terminal* if no later CNOT touches either endpoint. If a
degree-two wire is an endpoint of such an edge, choose an input Pauli
axis that commutes through its first incident CNOT. Its image is then
supported only on the terminal pair and is a one-CNOT image. The exact
Gram classification in global_lower_bound5_spectral_gram.md proves:

- A terminal edge on physically adjacent wires cannot have a degree-two
  endpoint.
- A terminal edge on opposite wires cannot have **both** endpoints of
  degree two. Each selected input axis would be off-diagonal and would
  commute with \(D^2\), forcing both \(V_4^2\) eigenspace dimensions to
  be divisible by four, contrary to 10 and 6.

The latter argument uses the same multiplicity step as the four-CNOT proof;
it remains valid even if other CNOTs occur earlier elsewhere.

## Exhaustive finite count

The script global_lower_bound6_terminal_filter.py enumerates all
\(6^5=7776\) chronological sequences of five *unoriented* wire pairs.
This covers all connectivities and allows both CNOT directions. It first
applies the degree condition, then the necessary full forward support of
\(U Z_jU^\dagger\) for each input wire \(j\), then the terminal-edge
obstructions. The exact output is
global_lower_bound6_terminal_filter_result.json (ledger run 190).

| Degree sequence | Degree feasible | Full forward support | Terminal excluded | Survive |
|---|---:|---:|---:|---:|
| \((4,2,2,2)\) | 720 | 216 | 176 | 40 |
| \((3,3,2,2)\) | 1860 | 792 | 488 | 304 |

These are counts of ordered pair schedules, not of locally inequivalent
circuits or verified diagonalizers. The surviving 344 schedules retain
arbitrary local gates and CNOT orientations. No assertion is made that any
survivor works. The next step needs an invariant for degree-two wires whose
second incident CNOT is followed by another CNOT on its partner; in that
case its selected observable can spread beyond the pair, so the terminal
Gram lemma no longer applies.
