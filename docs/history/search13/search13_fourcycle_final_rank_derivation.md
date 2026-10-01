# Four-cycle final-layer screen and fixed-target six-CNOT bound

Run 216 completes the degree-two interaction-graph screen left by run 213.
The construction fixes a five-CNOT parity prefix from the 196
one-mask-missing cases and the exact first `02/13` layer, including its
input local frames and pair-local endpoint map. The required final target
is the exact operator `F01 F23 C T` from
`search13_crosspair_second_layer_rank_derivation.md`, up to a right local
factor that leaves all operator Schmidt ranks unchanged.

The run enumerates the three undirected four-cycles on four wires. Each
cycle has 24 chronological edge orders and 16 direction choices, or 384
directed schedules; 1,152 schedules are covered in total. Arbitrary
one-qubit gates may be inserted around their four CNOTs.

Exact arithmetic in `Q(zeta_96)` gives the same operator Schmidt ranks
for all 51 distinct endpoint/mask/phase targets:

| Cut | Exact rank | Minimum crossing CNOTs |
|---|---:|---:|
| `01|23` | 5 | 3 |
| `02|13` | 16 | 4 |
| `03|12` | 16 | 4 |

The minimum crossing count is `ceil(log2(rank))`, because each CNOT
crossing a cut has rank two and operator Schmidt rank is submultiplicative.
Each two-qubit CNOT edge crosses exactly two of these three balanced cuts.
The three minimum counts sum to `3+4+4=11`, whereas four CNOTs contribute
only eight crossings in total. Thus **no four-CNOT final suffix topology**
can realize these fixed targets. In fact, any such suffix needs at least
`ceil(11/2)=6` CNOTs. The four-cycle schedules fail as a special case;
there was no compatible case for a numerical local-gate fit and no
13-CNOT candidate.

The target ranks on every one-wire cut are also four. They yield the
degree-two condition from run 213, but the balanced-cut counting argument
is stronger and covers all four-CNOT suffix graphs without enumerating
them.

This six-CNOT suffix bound applies only to the fixed exact first layer and
catalogued one-mask-missing phase-polynomial prefixes. It does not imply a
global 15-CNOT bound or exclude a 13-CNOT diagonalizer with a changed
first-layer operator, phase gauge, or output eigenspace basis.

A distinct next hypothesis is to change the first `02/13` layer itself so
that the transported omitted phase lowers at least one of the two rank-16
balanced cuts. An exact rank calculation for a parameterized first-layer
Cartan family should precede a bounded numerical joint fit.

Reproduce with:

```bash
python3 search13_fourcycle_final_rank.py
python3 experiment_log.py show 216
```
