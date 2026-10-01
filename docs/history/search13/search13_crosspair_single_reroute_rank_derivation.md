# One-CNOT reroutes of the final four-CNOT layer

The fixed-target construction starts with a five-CNOT parity prefix from
the 196 one-mask-missing cases, followed by the exact first `02/13` layer
with its pair-local endpoint map. The resulting target for the final
four-CNOT layer is the exact operator `F01 F23 C T`, as defined in
`search13_crosspair_second_layer_rank_derivation.md`. Right multiplication
by `C^dagger` is local across every cut and leaves the ranks unchanged.

Run 213 tests 16 directed final-layer CNOT schedules. Each replaces
exactly one of the chronological edges `01,01,23,23` with either direction
of a `03` or `12` CNOT. Arbitrary one-qubit gates are allowed between and
around these four CNOTs. The schedule still costs four CNOTs, giving 13
total with the five-CNOT prefix and four-CNOT first layer.

For a bipartition, a CNOT crossing the cut has operator Schmidt rank two;
a CNOT staying within one side and every one-qubit gate has rank one.
Therefore a schedule with `k` crossing CNOTs can implement only operators
of rank at most `2^k` across that cut. Run 213 calculates target ranks
exactly over `Q(zeta_96)` for all 51 distinct endpoint/mask/phase targets:

| Single-wire cut | Target rank in all 51 cases |
|---|---:|
| `0|123` | 4 |
| `1|023` | 4 |
| `2|013` | 4 |
| `3|012` | 4 |

Every tested one-edge reroute leaves one physical wire incident to only
one final-layer CNOT. Its rank cap is two, contradicting the target rank
four. Thus none of the `51 x 16` case/topology combinations can realize
the fixed target. No numerical local-rotation fit was run and no 13-CNOT
candidate emerged.

The same rank argument gives a useful necessary degree condition: any
four-CNOT final layer for these fixed targets must touch every wire at
least twice. Four CNOTs provide exactly eight incidences, so every wire
must have degree exactly two. Among four-edge multigraphs on four wires,
the remaining interaction shapes are a four-cycle or two doubled
disjoint edges. The disjoint `03/12` doublets were excluded by run 210;
the four-cycle is a distinct next topology to screen across balanced cuts
and then fit if it survives.

This conclusion fixes the exact first layer and target up to the
catalogued invariant input phase gauges. It is not a general lower bound
for all 13-CNOT cycle diagonalizers or for globally refitted eigenbases.

Reproduce with:

```bash
python3 search13_crosspair_single_reroute_rank.py
python3 experiment_log.py show 213
```
