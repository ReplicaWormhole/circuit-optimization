# Fixed middle-block three-CNOT obstruction

Let `F_ab(a,b) = exp(i a X_a X_b + i b Y_a Y_b)`. In the exact 14-CNOT
diagonalizer, the first blocks `F02(pi/4,pi/8)` and `F13(pi/4,pi/8)` are
disjoint and commute. The second blocks `F01(pi/8,pi/8)` and
`F23(pi/8,pi/8)` are also disjoint. Write the center local layer as
`L0 L1 L2 L3`. Since `L2` commutes with `F13` and `F01`, the full operator
can be regrouped, in matrix order, as

`U14 = F23 L2 [F01 (L0 L1 L3) F13] F02 L_input P6`.

The bracketed `M` acts only on wires `(0,1,3)` and the other pieces require
`6+2+2=10` CNOTs. A three-CNOT implementation of this **fixed exact**
`M` would therefore produce a 13-CNOT diagonalizer. An independent numerical
reconstruction of the regrouped four-qubit matrix agrees with the existing
exact 14-CNOT gate list to maximum entry error `5.92e-16`, up to global
phase. This numerical comparison checks the bookkeeping; the exact
matchgate identity is already established by the source construction.

Run 298 constructs `M` over `Q(zeta_96)`. For every ordered three-CNOT
schedule on the three physical pairs, it compares exact modular lower
bounds on mixed input/output realignment rank to a tensor-network mincut
upper bound. A CNOT factors as

`|0><0| tensor I + |1><1| tensor X`,

so its control and target factors have a dimension-two internal bond. Wire
segments also have dimension two, and arbitrary one-qubit gates can be
absorbed into adjacent factors. Thus a cut of `c` bonds caps realignment
rank at `2^c`, regardless of CNOT direction.

All `3^3=27` ordered pair schedules are excluded. The singleton-wire
operator-Schmidt ranks of `M` are all four, forcing any three-CNOT schedule
to use the triangle (each pair once). For each of its six orders, one of
mixed masks `13`, `19`, or `25` has target rank eight modulo both 97 and
193, hence rank eight over `Q(zeta_96)`, while the schedule's factor graph
has a two-bond cut and rank cap four. Independent NetworkX cuts reproduced
every recorded violation, and an independent modular elimination reproduced
all 31 target rank entries.

This excludes a three-CNOT synthesis of the specified `M` with arbitrary
one-qubit gates. It does not exclude other regroupings, intermediate/output
gauges, or a 13-CNOT diagonalizer with a different architecture.

Reproduce with `python3 search13_middle_threequbit.py` and inspect
`search13_middle_threequbit_result.json` (ledger run 298; supersedes the
pre-computation import failure in run 296).
