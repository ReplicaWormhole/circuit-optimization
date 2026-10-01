# Exact screen of a disjoint `03/12` second-layer reroute

This experiment keeps the exact first `02/13` layer (including its input
local frames and the pair-local endpoint permutation) and the exact tail
target up to the catalogued invariant input phase gauges. It replaces the
final disjoint `01/23` pair layer
with a disjoint `03/12` pair layer, allowing arbitrary local gates and
Cartan angles in the replacement blocks. The five-CNOT prefixes are the
196 one-mask-missing cases from runs 198 and 200.

Write the exact first-layer operator after the five-CNOT prefix as `A` and
the omitted parity rotation as `D`. The phase transported through that
layer is `T=A D A^dagger`. If `C` is the exact center local layer and
`F=F01 F23` is the original final disjoint pair layer, the full operator
that the replacement second layer must realize is

`S = F C T C^dagger`.

The right factor `C^dagger` is local across `03|12`, so it does not change
operator Schmidt rank across that cut. The exact screen tests `F C T`.
Any product of a `03` two-qubit block and a `12` two-qubit block has rank
one across `03|12`; this remains true with arbitrary Cartan angles and
one-qubit frames. Rank one is therefore necessary for this reroute.

Run 210 evaluates all 196 cases, deduplicated to 51 endpoint/mask/phase
operators and 20 transported generators. It uses exact arithmetic in
`Q(zeta_96)` and finds a nonzero `2x2` realignment minor across `03|12`
for **every** transported phase and **every** full operator `F C T`.
None has rank one. Independent indexing checks give rank one for identity
and `CX03 CX12`, and rank greater than one for `CX01`.

There is no compatible fixed-first-layer target for a numerical fit with
disjoint `03/12` second blocks, so this run produced no 13-CNOT candidate.
This does not exclude a globally refitted first layer, a mixed CNOT pattern
within the four-CNOT final layer, or a different eigenbasis inside the
degenerate shift eigenspaces.

A distinct next topology is to reroute only one of the four final-layer
CNOTs to `03` or `12`, preserving four final-layer CNOTs but allowing a
non-disjoint interaction graph. An exact cut-rank screen of those discrete
topologies should precede any local-rotation fit.

Reproduce with:

```bash
python3 search13_crosspair_second_layer_rank.py
python3 experiment_log.py show 210
```
