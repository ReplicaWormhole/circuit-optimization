# Finite same-eigenvalue row gauges below 14 CNOTs

The exact 14-CNOT diagonalizer has output eigenvalue powers
`(0,3,1,1,2,0,2,2,0,2,0,0,3,1,3,0)` in units of `i`. Any unitary
mixing rows with the same power preserves diagonalization. Run 219 tests
all 27 unordered equal-eigenvalue row pairs with one row swap or one real
Givens rotation by `pi/8`, `pi/4`, or `pi/2`: 108 finite gauges in all.

Two splits are examined. After the six-CNOT parity prefix, the full
eight-CNOT exact tail is `F01 F23 C F02 F13 L`; a 13-CNOT circuit with
the fixed prefix would require a seven-CNOT implementation of its gauged
version. After the first `F02/F13` layer, the final four-CNOT target is
`F01 F23 C`; a 13-CNOT circuit with the first ten CNOTs fixed would
require a three-CNOT implementation of its gauged version.

Run 219 computes seven operator Schmidt ranks exactly in `Q(zeta_96)` and
compares them with every undirected three- or seven-CNOT interaction
multigraph. None of the 108 gauged final-layer operators passes a
three-CNOT graph screen. In particular, each has rank at least 13 across
both `02|13` and `03|12`, while three CNOTs give rank at most eight across
any one cut. The fixed-first-layer, three-CNOT-final route is excluded for
this finite gauge family.

All 108 gauged full tails pass the necessary seven-CNOT graph screen,
with 108 compatible undirected multigraphs per gauge. This does not give
a synthesis. The row swap `(8,11)` yields the weakest full-tail balanced
rank profile in the scan, `(12,16,14)`. A `pi/4` Givens rotation on the
same rows gives `(16,13,16)`. Independent numerical reconstruction of the
rational gate list agrees with the block formula to `9.5e-16` up to global
phase; applying either gauge leaves the `V4` off-diagonal error at
`6.7e-16`.

Runs 220 and 222 make bounded seven-CNOT fixed-unitary fits on two graph
and edge-order choices closest to the eight-CNOT exact tail. They use
arbitrary SU(2) layers before, between, and after all CNOTs and maximize
phase-insensitive process overlap. The row swap has two seeds per topology;
the Givens rotation has one. Each fit is limited to 250 L-BFGS-B
iterations. The best checked 13-CNOT candidate has maximum off-diagonal
error about `0.634`, so none is a diagonalizer. These are restricted
numerical failures and do not exclude the other rank-compatible graphs,
orders, directions, or output eigenspace gauges.

A distinct next hypothesis is a joint eigenspace gauge acting on two
row pairs or a continuous three-dimensional block within one degenerate
sector. An exact rank objective could first seek a full-tail balanced
rank at most eight on one cut, reducing the required seven-CNOT graph
family before another bounded synthesis fit.

Reproduce with:

```bash
python3 search13_output_eigenrow_gauge_rank.py
python3 search13_eigenrow_gauge_seven_fit.py --maxiter 250
python3 search13_eigenrow_givens_seven_fit.py --maxiter 250
python3 experiment_log.py show 219
python3 experiment_log.py show 220
python3 experiment_log.py show 222
```
