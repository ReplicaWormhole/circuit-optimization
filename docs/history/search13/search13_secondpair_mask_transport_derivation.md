# One omitted parity phase transported to the second pair layer

The exact 14-CNOT circuit has a six-CNOT parity prefix, first disjoint pair
blocks on `02` and `13`, a center local layer, then disjoint pair blocks on
`01` and `23`. This test fixes the first layer's exact operator, including
its input local frames, but permits its Cartan class to change when a
pair-local prefix endpoint map is absorbed. It also permits arbitrary
two-qubit operations in the disjoint second pair blocks.

Runs 198 and 200 left 196 path/pair/model cases: each has a five-CNOT
prefix, exactly one omitted required parity mask, and an endpoint whose
linear mismatch from the exact six-CNOT prefix acts on `02` or `13` only.
For a case, let `P` be that pair-local endpoint permutation and let `L`
be the exact input local frame. On each first pair, the new block operator
is `A = F L P`, with `F = exp(i*pi*XX/4+i*pi*YY/8)`. A missing parity
rotation with generator `Z_m` must cross this first layer to be placed in
the second layer. Its transported generator is `Q = A Z_m A^dagger` (with
the two first-pair factors assembled on four wires). The center local
layer preserves its wire-support pattern and operator Schmidt rank across
`01|23`.

Run 204 computes the two first-pair factors exactly in `Q(zeta_96)`. It
checks 30 distinct factor/end-point combinations. In all 196 cases, each
factor remains supported on both wires of its first pair. Therefore no
transported generator is a local `Z` or a `ZZ` rotation on one second pair,
even allowing center local rotations of its Pauli axes.

Run 206 checks the **full phase operator**, accounting for the possibility
that an exponential could factor even when its generator has broad support.
The 196 cases reduce to 51 distinct endpoint/mask/angle triples. Their
omitted phase coefficients modulo 16 are
`{1,2,3,5,6,7,9,10,11,14,15}` in units of `pi/8`; none is a trivial or
`pi` phase. For each transported phase, exact realignment across `01|23`
has a nonzero `2x2` minor. Hence every phase has operator Schmidt rank at
least two across that cut. A product of any `01` block and any `23`
block has rank one, and multiplying it by center local gates or existing
second-layer pair blocks preserves the rank. Thus none of the 51 phases
can be absorbed solely in the second disjoint pair layer.

This excludes the fixed-first-layer, pair-local-endpoint, one-mask
relocation ansatz. It is not a lower bound on general 13-CNOT
diagonalizers: changing the first layer operator, inserting a cross-pair
interaction, or changing the eigenspace gauge can evade the rank-one
condition. No numerical fit was run because the exact necessary condition
fails for every case.

The next distinct construction test is to replace one second-layer edge
with a cross-pair edge (`03` or `12`) and screen the resulting seven/eight
CNOT tail for exact cut-rank compatibility before local fitting.

Reproduce with:

```bash
python3 search13_secondpair_mask_transport.py
python3 search13_secondpair_phase_rank.py
python3 experiment_log.py show 204
python3 experiment_log.py show 206
```
