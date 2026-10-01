# The other three overlapping matchgate pairs cannot save a CNOT

The exact 14-CNOT circuit has a six-CNOT prefix followed by first-layer
blocks `F02`, `F13`, center locals `L0 L1 L2 L3`, and final-layer blocks
`F01`, `F23`, where

`Fct(a,b) = exp(i a Xc Xt + i b Yc Yt)`.

The first layer's blocks commute with each other, as do the final layer's.
Choosing one block from each layer gives four possible overlapping three-wire
middle factors. Run 298 already excluded a three-CNOT synthesis of
`F01 L1 F13` (with its other center locals included). The other three,
after commuting center locals on nonshared wires to the exterior, are

| Factor | Physical wires | Three-CNOT schedules excluded |
|---|---|---:|
| `F01 L0 F02` | `0,1,2` | 27 of 27 |
| `F23 L2 F02` | `0,2,3` | 27 of 27 |
| `F23 L3 F13` | `1,2,3` | 27 of 27 |

Each of these fixed factors currently uses four CNOTs. If one admitted three,
the other first and final blocks could remain as two-CNOT factors and the
whole circuit would use `6+2+3+2=13` CNOTs. The exact screen excludes this
specific pairwise compression route.

## Exact screen

`search13_sibling_middle_blocks.py` constructs the three 8-by-8 target
matrices over `Q(zeta_96)`. For each of 31 complementary mixed input/output
cuts and each prime 97 and 193, it evaluates the target matrix at a primitive
96th root and computes its rank. Every denominator is invertible at both
primes. A nonzero finite-field minor gives a characteristic-zero rank lower
bound. Boundary axes are ordered `(out0,out1,out2,in0,in1,in2)` using each
factor's physical wire order.

The rank across each singleton physical wire is four (canonical masks
`9,45,27`). A three-CNOT circuit needs at least two CNOTs touching each
wire, so its three pair edges must form a triangle. For each of the six
triangle orders, a mixed cut has target rank eight but the two-dimensional
CNOT-factor tensor network has a two-bond cut, limiting circuit rank to
four. The script records a witness cut for each order and also excludes the
21 nontriangle schedules directly. Both CNOT directions and arbitrary
one-qubit gates are covered by the factor-network bound.

An independent floating-point construction from matrix exponentials agrees
with the three exact-field targets to maximum entry errors
`6.79e-16`, `6.47e-16`, and `6.13e-16`, respectively. Independent NetworkX
minimum cuts agree with all 837 brute-force schedule/cut capacities.
`search13_sibling_middle_blocks_result.json` contains every modular rank and
schedule witness. The ledger entry is run 311.

This proves a **fixed-operator, fixed-regrouping no-go**. It does not rule out
joint compression across more than two blocks, a changed intermediate or
output eigenbasis, a changed prefix, or a 13-CNOT diagonalizer elsewhere.
