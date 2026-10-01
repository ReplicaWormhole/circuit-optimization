# Exact output-row-gauge screens near the eight-CNOT tail

These runs target the exact 14-CNOT diagonalizer in
`topology14_exact_matchgate_rational.json`.  Its first six CNOTs are held
fixed.  The target of the seven-CNOT tail is $ST$, where $T$ is the
original exact eight-CNOT tail and $S$ permutes its output rows.  Any
such $S$ preserves diagonalization, with correspondingly permuted
eigenvalue labels.  The question here is whether a specified **pair
schedule** of seven CNOTs with arbitrary one-qubit gates could realize
that fixed target.  These screens do not search every output unitary
inside degenerate eigenspaces or every 13-CNOT architecture.

## Exact necessary condition

For an eight-bit mask $m$, regard the $16\times16$ target as a tensor
with four output and four input binary indices, realign by $m$, and
take its matrix rank.  There are 127 complementary mixed cuts.  The
entries of $T$ lie in $\mathbb Q(\zeta_{96})$.  Reducing its entries
modulo 97 or 193 at a primitive 96th root gives a rank **lower bound**
on the characteristic-zero rank.  A seven-CNOT circuit has a tensor
network with dimension-two wire bonds and one dimension-two factor bond
per CNOT; its realignment rank is at most $2^{\lambda_m}$, where
$\lambda_m$ is the bond-graph mincut for $m$.  Therefore a modular
rank greater than $2^{\lambda_m}$ excludes that pair schedule for
every CNOT direction and all one-qubit gates.

The near-tail family contains **124** distinct seven-pair schedules:
delete one of the eight CNOTs of $T$ and change at most one retained
unoriented pair.  The distance-two family changes at most two retained
pairs and contains **1,632** schedules, including the 124.

## Results

| Output-row gauge family | Gauges tested | Near-tail survivor pairs |
|---|---:|---:|
| One output CNOT | 12 | 0 |
| Two noncancelling output CNOTs | 132 ordered sequences | 0 |
| Three noncancelling output CNOTs | 1,452 sequences, 626 distinct maps | 0 |
| Every invertible binary linear map $\mathrm{GL}(4,2)$ | 20,160 distinct maps | **0** |
| One arbitrary output-row transposition | 120 | 0 |
| Two disjoint output-row transpositions | 5,460 | 0 |
| One oriented output-row three-cycle | 1,120 | 0 |

The exhaustive linear-map test (run 315) uses only reduction modulo 97:
758,232 adaptive cut-rank calculations suffice to reject all
$20{,}160\times124=2{,}499{,}840$ gauge/schedule pairs.  Its breadth-first
enumeration starts with the identity row map and composes all 12 directed
output CNOT transvections until it has exactly 20,160 maps.  An affine
translation of the four output bits is an output local $X$ layer, so
it leaves every mixed-cut rank unchanged and does not alter this
particular rank screen.

For distance at most two (run 316), all 20,160 linear maps leave **6,096**
gauge/schedule pairs across **1,368** gauges.  Every survivor has
distance exactly two.  The surviving schedules comprise only **13**
distinct unoriented pair orders; four already survive for the identity
gauge.  The additional nine were numerically fitted in run 319 with
their shortest compatible linear gauges.  All checked 13-CNOT circuits
failed the independent cycle diagonalization test; the best maximum
off-diagonal error was about 0.459.  Run 321 additionally fitted 18
topologically diverse identity-gauge survivors at distances three to
seven from a deletion; none passed the independent checker, and the
best maximum off-diagonal error was about 0.426.  These bounded
numerical failures do not exclude those schedules.

Run 322 extends the near-tail screen to all 1,120 oriented three-cycles
of output rows and rejects every one of their 124 schedule choices using
modulo-97 mixed-cut ranks.  Run 323 broadens the same gauges to the
1,632 distance-at-most-two schedules.  It leaves 416 gauge/schedule
pairs over 104 gauges and only four distinct schedules, indices 506,
530, 922, and 943 in its result file.  These four had already survived
with the identity output gauge.  Thus these particular three-cycle
gauges reveal no new pair order in this distance-two neighborhood.
Run 326 selected one of those three-cycles by the smallest sum of its
127 modular cut ranks and another with all three rows in one
eigenvalue sector.  Process-overlap and label-free fits on the two
representative surviving pair orders, each from deletion-warm and random
starts, failed the independent checker.  The best maximum
off-diagonal error was about 0.992.  Rank-profile selection is a
heuristic, and this numerical result excludes no surviving schedule.

Scripts/results: `search13_output_cx_near_tail_rank.py`,
`search13_two_output_cx_near_tail_rank.py`,
`search13_three_output_cx_near_tail_rank.py`,
`search13_all_linear_output_near_tail_rank.py`,
`search13_all_linear_output_distance2.py`,
`search13_row_swap_near_tail_rank.py`,
`search13_double_row_swap_near_tail_rank.py`,
`search13_threecycle_near_tail_rank.py`,
`search13_threecycle_distance2_rank.py`,
`search13_threecycle_survivor_fit.py`, and
`search13_linear_gauge_nine_fit.py`, with matching `_result.json` files.

The established global interval remains $6\leq C_{\min}\leq14$.
The exact exclusions above concern the fixed six-CNOT prefix, the
specified output-row gauge families, and the specified near-tail
seven-CNOT pair schedules only.
